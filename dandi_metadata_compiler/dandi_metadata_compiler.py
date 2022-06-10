from pathlib import Path
from argparse import ArgumentParser
import datetime
import json
import os

from tqdm import tqdm
from zetastitcher import FileMatrix, InputFile
from pyometiff import OMETIFFWriter
import yaml


def main():
    parser = ArgumentParser(description='Create DANDI XML from tiff files')
    parser.add_argument("input_dir", metavar="input_dir", type=str, nargs=1, help="Input directory")
    parser.add_argument("-c", "--config", metavar="config.yml", type=str, nargs=1, help="Config file, required",
                        required=True, dest="config")
    parser.add_argument("-o", "--output", metavar="output_dir", type=str, nargs=1,
                        help="Output directory, defaults to input directory",
                        required=False, default=None, dest="output_dir")
    parser.add_argument("--noxml", action="store_true", help="Do not create xml file", required=False)
    parser.add_argument("--nojson", action="store_true", help="Do not create json file", required=False)
    parser.add_argument("--nosymlinks", action="store_true", help="Do not create symlinks", required=False)
    parser.add_argument("--nochunks", action="store_true", default=False, help="dot not process chunks", required=False)
    parser.add_argument("--nomips", action="store_true", default=False, help="dot not process mips", required=False)
    parser.add_argument("--nofused", action="store_true", default=False, help="dot not process fused", required=False)

    args = parser.parse_args()

    input_dir = Path(args.input_dir[0])
    config_file = Path(args.config[0])
    output_dir = Path(args.output_dir[0]) if args.output_dir is not None else input_dir

    write_xml = not args.noxml
    write_json = not args.nojson
    make_symlinks = not args.nosymlinks

    process_chunks = not args.nochunks
    process_mips = not args.nomips
    process_fused = not args.nofused

    if not input_dir.exists():
        print("Input dir does not exist")
        exit(1)
    if not config_file.exists():
        print("Config file does not exist")
        exit(1)

    metacomp = DandiMetadataCompiler(input_dir=input_dir,
                                     config_file=config_file,
                                     output_dir=output_dir,
                                     process_chunks=process_chunks,
                                     process_mips=process_mips,
                                     process_fused=process_fused,
                                     write_xml=write_xml,
                                     write_json=write_json,
                                     make_symlinks=make_symlinks)
    metacomp.process_dir()


class DandiMetadataCompiler:
    def __init__(self,
                 input_dir: Path,
                 config_file: Path,
                 output_dir: Path = None,
                 output_json: bool = False,
                 process_chunks: bool = True,
                 process_mips: bool = True,
                 process_fused: bool = True,
                 write_xml: bool = True,
                 write_json: bool = True,
                 make_symlinks: bool = False,
                 skip_existing: bool = False):

        self.input_dir = input_dir
        self.path_parse_dict = self._parse_path()

        self.config_file = config_file
        self.output_json = output_json

        self.skip_existing = skip_existing
        self.stitchfile_path = self.input_dir.joinpath("stitch.yml")
        self.filematrix = FileMatrix(str(self.stitchfile_path))
        self.config_dict = self._parse_config_file(self.config_file)

        if output_dir is None:
            self.output_dir = self.input_dir
        else:
            self.output_dir = output_dir

        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True)

        self.write_xml = write_xml
        self.write_json = write_json
        self.make_symlinks = make_symlinks

        self.process_chunks = process_chunks
        self.process_mips = process_mips
        self.process_fused = process_fused

    @staticmethod
    def _get_value(key,
                   d: dict,
                   default_value=None):
        """
        Get value from a dict, if key is not in the dict, return default_value
        """
        if key not in d:
            return default_value
        else:
            return d[key]

    def _parse_config_file(self, config_file_fpath: Path) -> dict:
        yml_dict = self._read_yml(config_file_fpath)
        config_dict = {
            "Modality": self._get_value("Modality", yml_dict, "SPIM"),
            "PixelSize": self._get_value("PixelSize", yml_dict, [3.6, 3.6, 3.6]),
            "PixelSizeUnit": str(self._get_value("PixelSizeUnit", yml_dict, "µm")),
            "BitDepth": int(self._get_value("BitDepth", yml_dict, 16)),
            "Subject": str(self._get_value("Subject", yml_dict, self.path_parse_dict["subject"])),
            "Species": str(self._get_value("Species", yml_dict, "Homo Sapiens")),
            "BodyPart": str(self._get_value("BodyPart", yml_dict, "BRAIN")),
            "BodyPartDetails": str(self._get_value("BodyPartDetails", yml_dict, "")),
            "BodyPartDetailsExtended": str(self._get_value("BodyPartDetailsExtended", yml_dict, "")),
            "Environment": str(self._get_value("Environment", yml_dict, "exvivo")),
            # channel config
            "ExcitationWavelength": self._get_value("ExcitationWavelength", yml_dict, None),
            "EmissionWavelength": self._get_value("EmissionWavelength", yml_dict, None),
            "PhysicalUnit": str(self._get_value("PhysicalUnit", yml_dict, "nm")),
            "Name": str(self._get_value("Name", yml_dict, "")),
            "Fluor": str(self._get_value("Fluor", yml_dict, "")),
            # json common fields
            "ChunkTransformationMatrixAxis": self._get_value("ChunkTransformMatrixAxis", yml_dict, ["X", "Y", "Z"]),
            "BodyPartDetailsOntology": self._get_value("BodyPartDetailsOntology", yml_dict,
                                                       "http://purl.org/sig/ont/fma/fma68641"),
            "InstitutionAddress": self._get_value("InstitutionAddress", yml_dict,
                                                  "Via Nello Carrara 1, 50019 Sesto Fiorentino (FI), Italy"),
            "InstitutionName": self._get_value("InstitutionName", yml_dict,
                                               "LENS - European Laboratory for NonLinear Spectroscopy - University of Florence"),
            "InstitutionalDepartmentName": self._get_value("InstitutionalDepartmentName", yml_dict, "BioPhotonics Lab"),
            "Pathology": self._get_value("Pathology", yml_dict, "no pathology"),
            "SampleExtractionInstitution": self._get_value("SampleExtractionInstitution", yml_dict,
                                                           "Massachussetts General Hospital"),
            "SampleExtractionProtocol": self._get_value("SampleExtractionProtocol", yml_dict, "autopsy"),
            "SampleFixation": self._get_value("SampleFixation", yml_dict, None),
            "Sample_idx": int(self._get_value("Sample_idx", yml_dict, self.path_parse_dict["sample_idx"])),
            "AcquisitionDate": datetime.datetime.strptime(
                self._get_value("AcquisitionDate", yml_dict, self.path_parse_dict["date"]), "%Y%m%d").isoformat(),
            "Channels": yml_dict["Channels"],
            # "bftools_path": Path(yml_dict["paths_cfg"]["bftools_path"]),
            # "json_common": yml_dict["json_common"],
        }

        if config_dict["ExcitationWavelength"] is not None:
            config_dict["ExcitationWavelength"] = int(config_dict["ExcitationWavelength"])
        if config_dict["EmissionWavelength"] is not None:
            config_dict["EmissionWavelength"] = int(config_dict["EmissionWavelength"])

        if config_dict["Name"] == "":
            # if Name is not specified in config file, use the channel name obtained from the path
            # and get the wavelength from the channel name
            channel_wavelenght = self.path_parse_dict["channel"]
            channel_dict = config_dict["Channels"][channel_wavelenght]
            config_dict["Name"] = channel_dict["Name"]
            config_dict["Fluor"] = channel_dict["Fluor"]
            config_dict["ExcitationWavelength"] = channel_dict["ExcitationWavelength"]
            config_dict["EmissionWavelength"] = channel_dict["EmissionWavelength"]
            config_dict["PhysicalUnit"] = channel_dict["PhysicalUnit"]
        return config_dict

    @staticmethod
    def _read_yml(yml_fpath: Path) -> dict:
        """read a yml file, return a dict"""
        with yml_fpath.open(mode="r") as rfile:
            yml_dict = yaml.safe_load(rfile)
        return yml_dict

    @staticmethod
    def _get_stack_shape(stack_fpath: Path):
        return InputFile(stack_fpath).shape

    @staticmethod
    def _dump_dict_to_json(json_fpath: Path, json_dict: dict):
        if json_fpath.exists():
            json_fpath.unlink()
        with json_fpath.open(mode="w") as wfile:
            json.dump(json_dict, wfile, indent=4, sort_keys=True)

    @staticmethod
    def _make_symlink(in_fpath: Path,
                      link_fpath: Path):
        if link_fpath.is_file():
            assert link_fpath.is_symlink(), f"{str(link_fpath)} is not a symlink"
            os.remove(str(link_fpath))
        os.symlink(str(in_fpath), str(link_fpath))

    def _parse_path(self) -> dict:

        parent_path = self.input_dir.parent
        date_str, subject, sample_str, _, left_channel_str, _, right_channel_str = parent_path.name.split("_")
        sample_idx = int(sample_str)
        camera = self.input_dir.name.split("_")[1]
        if camera == "left":
            channel_str = left_channel_str
        elif camera == "right":
            channel_str = right_channel_str
        else:
            raise ValueError("Camera name is not left or right")
        channel_wavelength = int(channel_str)

        path_parse_dict = {
            "date": date_str,
            "subject": subject,
            "sample_idx": sample_idx,
            "channel": channel_wavelength,
            "camera": camera,
        }
        return path_parse_dict

    @staticmethod
    def _get_chunk_idx(self, chunk_fpath: Path) -> int:
        """get the chunk idx from the chunk fpath"""
        chunk_tiff_paths = list(chunk_fpath.parent.glob("x_*y_*z_.ome.tif"))
        # find the index of chunk name in the list of chunk tiff paths
        try:
            chunk_idx = chunk_tiff_paths.index(chunk_fpath)
        except ValueError:
            raise ValueError(f"{chunk_fpath} is not in {chunk_tiff_paths}")
        return chunk_idx

    def process_dir(self):
        """process input directory"""
        chunk_tiff_paths = sorted(list(self.input_dir.glob("x_*y_*z_*.ome.tif")))
        fused_path = self.input_dir.joinpath("fused.ome.tif")
        mip_path = self.input_dir.joinpath("mip.ome.tif")

        sub = self.config_dict["Subject"]
        ses = self.config_dict["Modality"]
        partdetails = self.config_dict["BodyPartDetails"]
        sample_idx = self.config_dict["Sample_idx"]
        sample = f"{partdetails}S{sample_idx:02d}"
        stain = self.config_dict["Name"]

        json_dict = {
            "PixelSize": self.config_dict["PixelSize"],
            "PixelSizeUnits": "um",
            "BodyPart": self.config_dict["BodyPart"],
            "BodyPartDetails": self.config_dict["BodyPartDetailsExtended"],
            "BodyPartDetailsOntology": self.config_dict["BodyPartDetailsOntology"],
            "Pathology": self.config_dict["Pathology"],
            "Environment": self.config_dict["Environment"],
            "SampleStaining": self.config_dict["Name"],
            "SampleExtractionProtocol": self.config_dict["SampleExtractionProtocol"],
            "SampleFixation": self.config_dict["SampleFixation"],
            "ChunkTransformationMatrixAxis": self.config_dict["ChunkTransformationMatrixAxis"],
            "InstitutionName": self.config_dict["InstitutionName"],
            "InstitutionAddress": self.config_dict["InstitutionAddress"],
            "InstitutionalDepartmentName": self.config_dict["InstitutionalDepartmentName"],
            "SampleExtractionInstitution": self.config_dict["SampleExtractionInstitution"],
        }

        ome_config_dict = {
            "AcquisitionDate": self.config_dict["AcquisitionDate"],
            "PhysicalSizeX": self.config_dict["PixelSize"][2],
            "PhysicalSizeXUnit": self.config_dict["PixelSizeUnit"],
            "PhysicalSizeY": self.config_dict["PixelSize"][1],
            "PhysicalSizeYUnit": self.config_dict["PixelSizeUnit"],
            "PhysicalSizeZ": self.config_dict["PixelSize"][0],
            "PhysicalSizeZUnit": self.config_dict["PixelSizeUnit"],
            "Channels": {
                self.config_dict["Name"]: {
                    "Name": self.config_dict["Name"],
                    "SamplesPerPixel": 1,
                    "ExcitationWavelength": self.config_dict["ExcitationWavelength"],
                    "ExcitationWavelengthUnit": self.config_dict["PhysicalUnit"],
                    "EmissionWavelength": self.config_dict["EmissionWavelength"],
                    "EmissionWavelengthUnit": self.config_dict["PhysicalUnit"],
                    "Fluor": self.config_dict["Fluor"],
                }
            }
        }

        # FUSED
        if self.process_fused and fused_path.is_file():
            fused_namestring = f"sub-{sub}_ses-{ses}_sample-{sample}_stain-{stain}_{ses}_fused"
            fused_xml_path = self.output_dir.joinpath(f"{fused_namestring}.xml")
            fused_symlink_path = self.output_dir.joinpath(f"{fused_namestring}.ome.tif")
            fused_ome_dict = ome_config_dict.copy()
            fused_ome_dict["Name"] = fused_namestring

            if self.write_xml:
                fused_array_shape = self._get_stack_shape(fused_path)
                fused_ome_writer = OMETIFFWriter(fpath=self.output_dir,
                                                 dimension_order="ZYX",
                                                 array=None,
                                                 metadata=fused_ome_dict,
                                                 arr_shape=list(fused_array_shape))
                fused_ome_writer.write_xml(fused_xml_path)
            if self.make_symlinks:
                self._make_symlink(in_fpath=fused_path, link_fpath=fused_symlink_path)

        # MIP
        if self.process_mips and mip_path.is_file():
            mip_namestring = f"sub-{sub}_ses-{ses}_sample-{sample}_stain-{stain}_{ses}_mip"
            mip_xml_path = self.output_dir.joinpath(f"{mip_namestring}.xml")
            mip_symlink_path = self.output_dir.joinpath(f"{mip_namestring}.ome.tif")
            mip_ome_dict = ome_config_dict.copy()
            mip_ome_dict["Name"] = mip_namestring
            if self.write_xml:
                mip_array_shape = self._get_stack_shape(mip_path)
                mip_ome_writer = OMETIFFWriter(fpath=self.output_dir,
                                               dimension_order="ZYX",
                                               array=None,
                                               metadata=mip_ome_dict,
                                               arr_shape=list(mip_array_shape))
                mip_ome_writer.write_xml(mip_xml_path)
            if self.make_symlinks:
                self._make_symlink(in_fpath=mip_path, link_fpath=mip_symlink_path)

        if self.process_chunks:
            for chunk_idx, chunk_tiff_path in enumerate(tqdm(chunk_tiff_paths)):
                chunk_namestring = f"sub-{sub}_ses-{ses}_sample-{sample}_stain-{stain}_chunk-{chunk_idx:02d}_{ses}"
                chunk_xml_out_path = self.output_dir.joinpath(f"{chunk_namestring}.xml")
                chunk_json_out_path = self.output_dir.joinpath(f"{chunk_namestring}.json")
                chunk_symlink_path = self.output_dir.joinpath(f"{chunk_namestring}.ome.tif")

                if self.skip_existing:
                    if chunk_xml_out_path.exists():
                        continue

                if self.write_xml:
                    chunk_ome_dict = ome_config_dict.copy()
                    chunk_ome_dict["Name"] = chunk_namestring

                    chunk_array_shape = self._get_stack_shape(chunk_tiff_path)
                    chunk_ome_writer = OMETIFFWriter(fpath=self.output_dir,
                                                     dimension_order="ZYX",
                                                     array=None,
                                                     metadata=chunk_ome_dict,
                                                     arr_shape=list(chunk_array_shape).copy())
                    chunk_ome_writer.write_xml(xml_fpath=chunk_xml_out_path)

                if self.write_json:
                    # write sidecar json
                    chunk_transform_matrix = self._get_chunktransformationmatrix(chunk_tiff_path)
                    chunk_json_dict = json_dict.copy()
                    chunk_json_dict.update({
                        "ChunkTransformationMatrix": chunk_transform_matrix,
                    })
                    self._dump_dict_to_json(chunk_json_out_path, chunk_json_dict)

                if self.make_symlinks:
                    self._make_symlink(in_fpath=chunk_tiff_path,
                                       link_fpath=chunk_symlink_path)

    def _get_chunktransformationmatrix(self, chunk_fpath: Path) -> list:
        """get the chunk transformation matrix from the chunk fpath"""
        chunk_key = "./" + chunk_fpath.name
        return [
            [0.0, 0.0, 0.0, int(self.filematrix.data_frame.loc[chunk_key, "Xs"])],
            [0.0, 0.0, 0.0, int(self.filematrix.data_frame.loc[chunk_key, "Ys"])],
            [0.0, 0.0, 0.0, int(self.filematrix.data_frame.loc[chunk_key, "Zs"])],
            [0.0, 0.0, 0.0, 1.0],
        ]


if __name__ == "__main__":
    main()
