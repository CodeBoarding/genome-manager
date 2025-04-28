#!/usr/bin/env python3

import os
import sys
from shlex import quote
import fileinput
from pathlib import Path
import argparse
import logging
import grp
import json
from platform import python_version
import shutil
import yaml
import glob
from pydantic import BaseModel, validator, FilePath, DirectoryPath, ValidationError
from typing import Optional, Union
from collections import defaultdict
import hashlib
import urllib.request
from urllib.error import HTTPError
from getpass import getuser
from generate_gtf_entry import YamlGeneCollection, YamlGeneModel
import generate_gtf_entry

__author__ = 'Rob Moccia'
__version__ = '0.2'

# global variables to define directory structure relative to top level
GENOMES_RELATIVE_PATH = Path('genomes')
USER_GENES_RELATIVE_PATH = Path('user_defined_genes')
GENOMES_CONFIG_DIR_RELATIVE_PATH = Path('.conf/genome-registry')
USER_GENES_CONFIG_DIR_RELATIVE_PATH = Path('.conf/user-registry')
MOUNTS_CONFIG_RELATIVE_PATH = Path('.conf/mounts.json')
LOG_DIR_RELATIVE_PATH = Path('.log')
TEMP_DIR_RELATIVE_PATH = Path('.tmp')
TEMP_DOWNLOAD_RELATIVE_PATH = Path(TEMP_DIR_RELATIVE_PATH, 'downloads')

# Custom exceptions
class DuplicateGenomeError(Exception):
    pass

class FileMismatchError(Exception):
    pass

class FileFormatError(Exception):
    pass

class InvalidSystemName(Exception):
    pass

class ChecksumMismatchError(Exception):
    pass

class MultipleMatchesError(Exception):
    pass

class NoMatchesError(Exception):
    pass

# Schemas
class GenomeFile(BaseModel):
    """Schema for a file that is part of a genome."""

    type: str
    default_system: str
    active_system: str=None
    path: dict[str, Path] # one each for e.g., hpc, aws, gpfs
    checksum: str=None
    source: Optional[str]
    content: Optional[str]
    parent: Optional[Path]
    
    @validator('type', pre=True)
    def validate_type(cls, value):
        valid_types = ['fasta', 'gtf', 'refflat', 'rrna_interval_list', 'gff3', 'yaml_gene_model']
        if value not in valid_types:
            raise ValueError(f'{value} is not a recognized type ({valid_types})')
        return value

    @validator('active_system', always=True)
    def set_active_system(cls, active_system: str, values: dict[str, any]) -> DirectoryPath:
        if active_system is None and 'default_system' in values:
            return values['default_system']
        else:
            return active_system
        
    @validator('checksum', always=True)
    def add_checksum(cls, val, values):
        try:
            if 'path' in values and 'active_system' in values:
                target = values['path'][values['active_system']]
                filesize = Path(target).stat().st_size
                if filesize > 100000:
                    return filesize
                with open(target, 'rb') as f:
                    file_hash = hashlib.md5()
                    while chunk := f.read(16384):
                        file_hash.update(chunk)
                    md5_hash = file_hash.hexdigest()
                try:
                    if val is not None:
                        if val != md5_hash:
                            raise ChecksumMismatchError(f'{target} checksum {md5_hash} does not match stored checksum {val}')
                        else:
                            logger.debug(f'{target} checksum {md5_hash} matches stored checksum {val}')
                except ChecksumMismatchError as e:
                    logger.exception(e)
                return md5_hash
        except Exception as e:
            logger.exception(f'values: {values}\n{e}')
            raise            

    @validator('source')
    def validate_source(cls, value):
        valid_sources = ['genome', 'transcriptome']
        if value is None:
            return
        if value not in valid_sources:
            raise ValueError(f'{value} is not a recognized source ({valid_sources})')
        return value

    # if the file is specified as type yaml_gene_model must validate that it
    # follows the specification for YamlGeneModel or YamlGeneCollection in generate_gtf_entry.py
    @validator('path', pre=True)
    def validate_path(cls, value, values):
        if values['type'] == 'yaml_gene_model' and 'active_system' in values:
            yaml_file = value[values['active_system']]
            try:
                logger.debug(f'trying to open {yaml_file}')
                with open(yaml_file, 'r') as f:
                        input_obj = yaml.load(f, Loader=yaml.CLoader)
                logger.debug(f'opened {yaml_file}')
            except Exception as e:
                logger.exception(
                    f'type = yaml_gene_model and {yaml_file} is not a valid YAML file\n{e}')
                raise
            try:
                if isinstance(input_obj, dict):
                    YamlGeneCollection.parse_obj(input_obj)
                elif isinstance(input_obj, list):
                    YamlGeneModel.parse_obj(input_obj[0])
                else:
                    raise ValidationError
            except NameError as e:
                logger.exception(f'No json config files could be found for {values["id"]}\n{e}')
                raise
            except ValidationError as e:
                logger.exception(
                    f'type = yaml_gene_model and {yaml_file} is not a valid YAML gene model specification\n{e}')
                raise
        return value

class GenomePath(BaseModel):
    """Schema for a directory that is part of a genome."""

    type: str
    default_system: str
    active_system: str=None
    # Path is a DirectoryPath, but validating as such will cause errors because not all mount points
    # are reachable on a given system. Existence is checked when adding via the register-genome process
    # so can use 'str' here safely.
    path: dict[str, Path] # one each for e.g., hpc, aws, gpfs
    source: Optional[str]
    content: Optional[str]
    parent: Optional[Path]
    
    @validator('type')
    def validate_type(cls, value):
        valid_types = ['star_index']
        if value not in valid_types:
            raise ValueError(f'{value} is not a recognized type ({valid_types})')
        return value

    @validator('active_system', always=True)
    def set_active_system(cls, active_system: str, values: dict[str, any]) -> DirectoryPath:
        if active_system is None and 'default_system' in values:
            return values['default_system']
        else:
            return active_system
        
    @validator('source')
    def validate_source(cls, value):
        valid_sources = ['genome', 'transcriptome']
        if value is None:
            return
        if value not in valid_sources:
            raise ValueError(f'{value} is not a recognized source ({valid_sources})')
        return value

class GenomeMetadata(BaseModel):
    """Schema for genome assembly metadata"""
    id: str
    species: str
    species_short: str
    release: int
    assembly: str
    assembly_type: str
    sequence_type: str

    @validator('assembly_type')
    def validate_assembly_type(cls, value):
        valid_assembly_types = ['pa', 'tl', 'primary_assembly', 'toplevel']
        if value not in valid_assembly_types:
            raise ValueError(f'{value} is not a valid assembly type ({valid_assembly_types})')
        return value

    @validator('sequence_type')
    def validate_sequence_type(cls, value):
        valid_sequence_types = ['dna', 'dna_rm', 'dna_sm']
        if value not in valid_sequence_types:
            raise ValueError(f'{value} is not a valid sequence type ({valid_sequence_types})')
        return value

class BaseGenome(BaseModel):
    """Schema for the core files representing a genome assembly (i.e., genome fasta and GTF)"""

    metadata: GenomeMetadata
    genome_fasta: GenomeFile
    gtf: GenomeFile
    description: Optional[str]

# Genome class attribute lists (useful for building/modifying a new object as a dictionary before passing to the schema)
GENOME_CLASS_MAIN_ATTR = ['transcriptome_fasta', 'refflat', 'rrna_interval_list', 'star_index']
GENOME_CLASS_BASE_ATTR = ['genome_fasta', 'gtf'] # located under 'base' key

class Genome(BaseModel):
    """Schema for representing a genome assembly, annotation, and associated files and metadata"""

    id: str
    default_system: str
    base: BaseGenome
    transcriptome_fasta: GenomeFile
    star_index: GenomePath
    refflat: GenomeFile
    rrna_interval_list: GenomeFile
    active_system: str=None
    description: Optional[str]

    def _attributes(self):
        """
        Helper function to return list of attributes with paths to be updated
        (e.g., when adding a new mountpoint)
        """
        return [self.base.genome_fasta, self.base.gtf, self.transcriptome_fasta, self.refflat,
                self.rrna_interval_list, self.star_index]

    def add_new_mountpoint(self, mountpoint: Union[str, bytes, os.PathLike], system_name: str,
                       verify: bool=True) -> None:
        """Add a new mountpoint to the path dictionary of all class files with that attribute"""
        logger.info(f'add_new_mountpoint called for {self.id}; mountpoint={mountpoint}; system_name={system_name}')
        for attribute in self._attributes():
            add_new_basepath(attribute, GENOMES_RELATIVE_PATH,
                             basepath=mountpoint, system_name=system_name,
                             verify=verify)
            logger.debug(f'added mountpoint for {self.id} {attribute.type}')

    # def remove_mountpoint(self, remove_system_name: str) -> None:
    #     """Remove a previously stored mountpoint under key 'remove_system_name'"""
    #     logger.info(f'remove_mountpoint called for {self.id}; system_name={remove_system_name}')
    #     found = False
    #     for attribute in self._attributes():
    #         if attribute.path.pop(remove_system_name, None):
    #             found = True
    #             logger.debug(f'removed mountpoint stored as system_name {remove_system_name} for {self.id} {attribute.type}')
    #             if attribute.active_system == remove_system_name:
    #                 attribute.active_system = None
    #                 logger.debug(f'set active_system to None for {self.id} {attribute.type}')
    #         else:
    #             logger.debug(f'system_name {remove_system_name} is not a key in {self.id} {attribute.type}')
    #     if found and self.active_system == remove_system_name:
    #         self.active_system = None
    #         logger.debug(f'set active_system to None for {self.id}')

    def propagate_active_system(self, active_system: str) -> None:
        """Push active system to the active_system attribute of all class attributes"""
        self.active_system = active_system
        self.base.genome_fasta.active_system = active_system
        self.base.gtf.active_system = active_system
        self.transcriptome_fasta.active_system = active_system
        self.star_index.active_system = active_system
        self.refflat.active_system = active_system
        self.rrna_interval_list.active_system = active_system

class GenomeCollection(BaseModel):
    """Schema for representing a collection of Genome objects"""

    genomes: dict[str, Genome]

    def get_genome_info(self):
        """
        Return a report of all genomes stored in the collection.
        """
        genome_info = defaultdict(dict)
        for id, genome in self.genomes.items():
            genome_info[genome.base.metadata.species] = {
                'id': id,
                'release': genome.base.metadata.release, 
                'assembly': genome.base.metadata.assembly}
        return genome_info

class UserDefinedGene(BaseModel):
    """Schema for representing a transcript(s) that can be added to a genome"""
    
    default_system: str
    active_system: str=None
    gene_model: dict[int, GenomeFile]
    fasta: GenomeFile
    id: str

    class Config:
        validate_assignment = True

    @validator('active_system', always=True)
    def set_active_system(cls, active_system: str, values: dict[str, any]) -> str:
        if active_system is None and 'default_system' in values:
            return values['default_system']
        else:
            return active_system
        
    @validator('gene_model')
    def validate_gene_model(cls, val: dict[int, GenomeFile], values):
        """Enforce the same gene_id in every YAML associated with the same UserDefinedGene"""
        if 'active_system' in values:
            ids = set()
            for genome_file in val.values():
                yaml_file = genome_file.path[values['active_system']]
                try:
                    yaml_obj = yaml.load(Path(yaml_file).open(), Loader=yaml.CLoader)
                except Exception as e:
                    logger.exception(f'failed to open {yaml_file}\n{e}')
                try:
                    current_id = YamlGeneModel.parse_obj(yaml_obj[0]).gene_id
                    ids.add(current_id)
                except NameError as e:
                    logger.exception(f'No json config files could be found for {genome_file}\n{e}')
                    raise
                try:
                    if len(ids) > 1:
                        raise ValueError(f'YAML files for different versions have different gene_id: {ids}')
                except Exception as e:
                    logging.exception(e)
                    raise
        return val

    @validator('fasta')
    def validate_fasta(cls, val: GenomeFile, values) -> GenomeFile:
        """Check that fasta file has only one entry named exactly the same as gene_id in the YAML model"""
        if 'gene_model' in values and 'active_system' in values:
            check_file = next(iter(values['gene_model'].values()))
            yaml_file = check_file.path[values['active_system']]
            yaml_obj = yaml.load(Path(yaml_file).open(), Loader=yaml.CLoader)
            gene_id = YamlGeneModel.parse_obj(yaml_obj[0]).gene_id
            filename = val.path[values['active_system']]
            description = []
            with open(filename, 'r') as f:
                for line in f:
                    if not line.startswith('>'):
                        continue
                    else:
                        description.append(line)
            if len(description) != 1:
                logger.error(
                    f'fasta file {filename} appears to have more than one sequence entry')
                raise ValueError('fasta file must have only 1 entry')
            else:
                fasta_id = description.pop().strip().lstrip('>')
                if fasta_id != gene_id:
                    logger.error(
                        f'sequence name in {filename} ({fasta_id}) does not match gene_id in YAML ({gene_id})')
                    raise ValueError('fasta file sequence name must match gene_id in YAML gene model')
            return val

    @validator('id')
    def validate_id(cls, val: str, values) -> str:
        if 'gene_model' in values and 'active_system' in values:
            latest_version = sorted(values['gene_model'].keys())[-1]
            check_file = values['gene_model'][latest_version]
            yaml_file = check_file.path[values['active_system']]
            yaml_obj = yaml.load(Path(yaml_file).open(), Loader=yaml.CLoader)
            if val != (gene_id := YamlGeneModel.parse_obj(yaml_obj[0]).gene_id):
                raise ValidationError(f'provided id ({val}) does not match gene_id in {yaml_file} ({gene_id})')
            return val

    def add_new_mountpoint(self, mountpoint: Union[str, bytes, os.PathLike], system_name: str,
                       verify: bool=True) -> None:
        """Add a new mountpoint to the path dictionary of all class files with that attribute"""
        logger.info(f'add_new_mountpoint called for {self.id}; mountpoint={mountpoint}; system_name={system_name}')
        add_new_basepath(self.fasta, USER_GENES_RELATIVE_PATH,
                         basepath=mountpoint, system_name=system_name,
                         verify=verify)
        logger.debug(f'added mountpoint for {self.id} fasta file')
        for version_number, model in self.gene_model.items():
            add_new_basepath(model, USER_GENES_RELATIVE_PATH,
                             basepath=mountpoint, system_name=system_name,
                             verify=verify)
            logger.debug(f'added mountpoint for {self.id} yaml file, version {version_number}')

    # def remove_mountpoint(self, remove_system_name: str) -> None:
    #     """Remove a previously stored mountpoint under key 'remove_system_name'"""
    #     logger.info(f'remove_mountpoint called for {self.id}; system_name={remove_system_name}')
    #     try:
    #         if self.default_system == remove_system_name:
    #             raise ValueError(f'Cannot remove default_system: {remove_system_name} is the default_system for {self.id}')
    #     except Exception as e:
    #         logger.exception(e)
    #         raise

    #     found = False
    #     if self.fasta.path.pop(remove_system_name, None):
    #         logger.debug(f'removed mountpoint stored as system_name {remove_system_name} for {self.id} fasta file')
    #         if self.fasta.active_system == remove_system_name:
    #             self.fasta.active_system = None
    #             logger.debug(f'set active_system to None for {self.id} fasta file')
    #     else:
    #         logger.debug(f'system_name {remove_system_name} is not a key in {self.id} fasta file')
    #     for version_number, model in self.gene_model.items():
    #         if model.path.pop(remove_system_name, None):
    #             found = True
    #             logger.debug(f'removed mountpoint stored as system_name {remove_system_name} for {self.id} yaml file, version {version_number}')
    #             if model.active_system == remove_system_name:
    #                 model.active_system = None
    #                 logger.debug(f'set active_system to None for {self.id} yaml file, version {version_number}')
    #         else:
    #             logger.debug(f'system_name {remove_system_name} is not a key in {self.id} yaml file, version {version_number}')
    #     if found:
    #         self.active_system = None
    #         logger.debug(f'set active_system to None for {self.id}')

    def propagate_active_system(self, active_system: str) -> None:
        """Push the active_system attribute to all class attributes"""
        self.active_system = active_system
        self.fasta.active_system = active_system
        for model in self.gene_model.values():
            model.active_system = active_system

    def get_version(self, version: int, system_name: str) -> Union[str, bytes, os.PathLike]:
        """
        Return path to a gene_model version by version number and system name.
        Latest version is represented by code -1.
        """
        version = int(version)
        if version < 0:
            version = sorted(self.gene_model.keys())[-1]
        return Path(self.gene_model[version].path[system_name])

    def add_version(self, yaml_file: Union[str, bytes, os.PathLike],
            system_name: str) -> tuple[int, Union[str, bytes, os.PathLike]]:
        """
        Add a new YAML gene_model to the gene_model dict, auto-incrementing the version
        Returns a tuple of the newly assigned version number as well as path where YAML
        was written which can be used to restore registry to pre-existing state on downstream
        failure.

        """
        latest_version = sorted(self.gene_model.keys())[-1]
        this_version = latest_version + 1
        # ensure consistent path for system_name by recovering it from a previous version
        try:
            prior_yaml_path = list(self.gene_model.values())[-1].path[system_name]
            registry_path = prior_yaml_path.parent
        except KeyError as e:
            logger.exception(f'{system_name} is not a valid system_name for {self.id}\n{e}')
            raise
        except Exception as e:
            logger.exception(f'failed to add new version to {self.id}\n{e}')
            raise
        yaml_dest = Path(registry_path,
            self.id + '_v' + str(this_version).zfill(2) + '.yaml')
        try:
            shutil.copy(yaml_file, yaml_dest)
            new_model = GenomeFile(default_system=self.default_system,
                active_system=system_name,
                path={system_name: yaml_dest}, type='yaml_gene_model')
            for version, next_model in self.gene_model.items():
                if new_model.checksum == next_model.checksum:
                    logger.error(f'checksum matches model stored in version {version}')
                    raise ValueError(f'YAML gene model is identical to a previously stored version')
            self.gene_model[latest_version + 1] = new_model
        except Exception as e:
            logger.exception(f'failed to add {yaml_file} as new version for {self.id}\n{e}')
            if yaml_dest.exists():
                yaml_dest.unlink()
                logger.debug(f'ERROR RECOVERY: undid addition of {yaml_dest} to registry')
            else:
                logger.debug(
                    'ERROR RECOVERY: no files were written so no action necessary to restore registry to pre-existing state')
            raise
        return latest_version + 1, yaml_dest

class UserDefinedGeneCollection(BaseModel):
    """Schema for representing all available modifications that can be added to a genome"""

    modifications: dict[str, UserDefinedGene]

class Mountpoints(BaseModel):
    """
    Schema for storing the mount points that have been added to a genome registry
    NOT CURRENTLY USED
    """

    default_system_name: str
    mounts: dict[str, Path]

    def __str__(self):
        column_width = 20
        res = [f'{"system_name": <{column_width}}{"mount point": <{column_width}}']
        for system_name, mount in self.mounts.items():
            res.append(f'{system_name: <{column_width}}{mount!s: <{column_width}}')
        return '\n'.join(res)
    
## Utility helper functions ##
def find_active_system(registry_path: Union[str, bytes, os.PathLike]) -> str:
    """
    Uses the mounts.config file to look up the system_name associated with a provided registry_path
    """
    system_found = False
    mount_config = load_mount_config(registry_path)
    for sysname, regpath in mount_config.mounts.items():
        if Path(regpath) == Path(registry_path):
            system_name = sysname
            system_found = True
            break
    if system_found:
        return system_name
    else:
        # this should not be possible, but include for debugging purposes just in case
        raise ValueError(f'could not find system_name for registry_path: {registry_path}')

def abbreviate_species(species: str) -> str:
    """
    Convert a species name in the form of <genus>_<species> to a 4-letter string consisting of
    the first letter of the genus and the first 3 letters of the species
    """
    fields = species.lower().split("_")
    return fields[0][:1] + fields[-1][:3]

def copy_with_logging(src: Union[str, bytes, os.PathLike], dest: Union[str, bytes, os.PathLike]
        ) ->  Union[str, bytes, os.PathLike]:
    """Wrapper around shutil.copy and shutil.copytree that adds logging"""
    src = Path(src)
    dest = Path(dest).resolve()
    if src.is_file():
        try:
            logger.info(f'copying {src} to {dest}')
            if not dest.parent.exists():
                dest.parent.mkdir(parents=True)
            if not dest.exists():
                write_path = shutil.copy(src, dest)
                logger.info(f'successfully copied {src} to {dest}')
            else:
                logger.error(f'attempted to overwrite existing file {dest} when copying {src}')
                raise FileExistsError(f'{dest} already exists')
        except Exception as e:
            logger.exception(f'failed to copy {src} to {dest}\n{e}')
            raise
    elif src.is_dir():
        try:
            logger.info(f'copying directory tree at {src} to {dest}')
            write_path = shutil.copytree(src, dest)
            logger.info(f'successfully copied directory tree at {src} to {dest}')
        except Exception as e:
            logger.exception(f'failed to copy directory tree at {src} to {dest}\n{e}')
            raise
    else:
        try:
            raise TypeError(f'{src} is not a file or a directory')
        except Exception as e:
            logger.exception(f'failed to copy {src} to {dest}\n{e}')

    return write_path

def globber(dir: Union[str, bytes, os.PathLike], pattern: str) -> Path:
    """Wrapper around glob that raises exceptions if no matches or more than one match"""
    os.chdir(dir)
    matches = glob.glob(pattern)
    if len(matches) == 1:
        return Path(dir, matches.pop())
    elif len(matches) > 1:
        raise MultipleMatchesError
    else:
        raise NoMatchesError

def glob_genome_files(dir: Union[str, bytes, os.PathLike]) -> dict[str, str]:
    """
    Uses glob.glob to find genome files in a given directory. This is a convenience function enabling
    new genomes to be built by providing all necessary genome files in a single directory rather than
    having to name them individually as arguments to register-genome.

    If new genome file types are added to the registry, they must also be added to GLOB_DICT so that they
    are searched for when running register-genome.
    """
    GLOB_DICT = {
        'genome_fasta': '*dna*.[fa][fasta].gz',
        'gtf': '*.gtf.gz',
        'transcriptome_fasta': '*.transcriptome*.[fa][fasta].gz',
        'refflat': '*.refflat',
        'rrna_interval_list': '*.rrna',
        'star_index': 'star-index*'
    }
    genome_files = {}
    for filetype, pattern in GLOB_DICT.items():
        try:
            res = globber(dir, pattern)
            genome_files[filetype] = res
        except Exception as e:
            logger.exception(f'required genome file type {filetype} not found in {dir}\n{e}')
            raise
    return genome_files

def format_assembly_name(assembly: str) -> str:
    """
    Helper function to reformat long or unwieldy genome assembly names into a shorter format
    suitable to be used as an ID
    """
    if assembly.lower().startswith("macaca"):
        fields = assembly.lower().split("_")
        assembly = fields[0][:1] + fields[1][:3] + fields[2]
    # else:
    #     assembly = assembly.replace(".", "")
    return assembly

def load_mount_config(registry_path: Union[str, bytes, os.PathLike]) -> Mountpoints:
    """Helper function to load mount points from config file"""
    mount_config_path = Path(registry_path, MOUNTS_CONFIG_RELATIVE_PATH).resolve()
    mounts = Mountpoints.parse_file(mount_config_path)
    logger.info(f'loaded mountpoint config file from {mount_config_path}')
    return mounts

def write_mount_config(registry_path: Union[str, bytes, os.PathLike],
                      mount_config: Mountpoints) -> None:
    """Helper function to save mount points to config file"""
    with open(os.path.join(registry_path, MOUNTS_CONFIG_RELATIVE_PATH), 'w') as f:
        f.write(mount_config.json())

def set_active_system_genome(genome_collection: dict, system_name: str) -> dict:
    for genome_name in genome_collection['genomes']:
        genome_collection['genomes'][genome_name]['active_system'] = system_name
        for key in GENOME_CLASS_BASE_ATTR:
            genome_collection['genomes'][genome_name]['base'][key]['active_system'] = system_name
        for key in GENOME_CLASS_MAIN_ATTR:
            genome_collection['genomes'][genome_name][key]['active_system'] = system_name
    return genome_collection

def load_genome(registry_file: Union[str, bytes, os.PathLike], system_name: str) -> GenomeCollection:
    """
    Parse a GenomeCollection object from a registry JSON filepath
    
    This function is required because the active_system needs to be set for the schema prior to validation
    or else it will attempt to load files on alternate mount points that aren't reachable from the system in
    use.
    """
    try:
        with open(registry_file, 'r') as f:
            genome_collection = json.load(f)
        genome_collection = set_active_system_genome(genome_collection=genome_collection, system_name=system_name)
        # for genome_name in genome_collection['genomes']:
        #     genome_collection['genomes'][genome_name]['active_system'] = system_name
        #     for key in GENOME_CLASS_BASE_ATTR:
        #         genome_collection['genomes'][genome_name]['base'][key]['active_system'] = system_name
        #     for key in GENOME_CLASS_MAIN_ATTR:
        #         genome_collection['genomes'][genome_name][key]['active_system'] = system_name
        model = GenomeCollection(**genome_collection)
        # model = GenomeCollection.parse_file(Path(registry_file))
        # for genome in model.genomes.values():
        #     genome.propagate_active_system(system_name)
        logger.info(f'loaded GenomeCollection from {Path(registry_file).resolve()}')
        return model
    except Exception as e:
        logger.exception(f'failed to load from {Path(registry_file).resolve()}\n{e}')
        raise

def write_genome_files(genome_dir: Union[str, bytes, os.PathLike], genome_fasta: Union[str, bytes, os.PathLike],
        gtf: Union[str, bytes, os.PathLike], transcriptome_fasta: Union[str, bytes, os.PathLike],
        refflat: Union[str, bytes, os.PathLike], rrna_interval_list: Union[str, bytes, os.PathLike],
        star_index: Union[str, bytes, os.PathLike]) -> dict[str, Union[str, bytes, os.PathLike]]:
    """
    Write genome files to the registry enforcing a consistent directory structure.
    Currently implemented specifically for Ensembl genomes.
    """
    genome_fasta_path = copy_with_logging(genome_fasta, Path(genome_dir, 'source', genome_fasta.name))
    gtf_path = copy_with_logging(gtf, Path(genome_dir, 'source', gtf.name))
    transcriptome_fasta_path = copy_with_logging(transcriptome_fasta, Path(genome_dir, 'derived', transcriptome_fasta.name))
    refflat_path = copy_with_logging(refflat, Path(genome_dir, 'derived', refflat.name))
    rrna_interval_list_path = copy_with_logging(rrna_interval_list, Path(genome_dir, 'derived', rrna_interval_list.name))
    star_index_path = copy_with_logging(star_index, Path(genome_dir, star_index.name))
    return {
        'genome_fasta': genome_fasta_path,
        'gtf': gtf_path,
        'transcriptome_fasta': transcriptome_fasta_path,
        'refflat': refflat_path,
        'rrna_interval_list': rrna_interval_list_path,
        'star_index': star_index_path
        }

def gene_model_from_yaml(yaml_file: Union[str, bytes, os.PathLike]) -> YamlGeneModel:
    """Create YamlGeneModel object from a YAML file"""
    try:
        yaml_obj = yaml.load(Path(yaml_file).open(), Loader=yaml.CLoader)
    except:
        raise
    if isinstance(yaml_obj, dict):
        raise ValueError(f'{yaml_file} has top-level key(s) -- consider register-gene-collection')
    elif not isinstance(yaml_obj, list):
        raise ValueError(f'{yaml_file} is not a valid YAML gene model')
    else:
        if len(yaml_obj) > 1:
            raise ValueError(f'{yaml_file} contains {len(yaml_obj)} entries -- remove extra entries or re-format as collection')
        model = YamlGeneModel.parse_obj(yaml_obj[0])
        return model

def parse_genome_metdata_file(json_args: Union[str, bytes, os.PathLike]) -> dict:
    """Read genome metdata from JSON file and return as dictionary with id converted to lowercase"""
    try:
        with open(json_args, 'r') as f:
            params = json.load(f)
        # convert id to lower by convention to enable case-insensitive searching
        params['id'] = params['id'].lower()
    except FileNotFoundError as e:
        logger.exception(f"failed to load {Path(json_args).resolve()}\n{e}")
        raise
    except KeyError as e:
        logger.exception(f'{Path(json_args).resolve()} is not formatted properly\n{e}')
        raise
    except AttributeError as e:
        logger.exception(f'failed to convert id to lowercase\n{e}')
        raise
    except Exception as e:
        logger.exception(e)
        raise
    return params

def humansize(nbytes: int) -> str:
    """
    Helper function to convert file length into human readable
    Modified from: https://stackoverflow.com/questions/14996453/python-libraries-to-calculate-human-readable-filesize-from-bytes
    as posted by user nneonneo
    """
    suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
    nbytes = int(nbytes)
    i = 0
    while nbytes >= 1024 and i < len(suffixes)-1:
        nbytes /= 1024.
        i += 1
    f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
    return '%s %s' % (f, suffixes[i])

def validate_user_gene_file(filename:  Union[str, bytes, os.PathLike]) -> None:
    """
    Validate that a file ends with a newline. This is necessary for files that will be
    concatenated prior to returning like those for user-defined genes.
    """
    with open(filename, 'r') as f:
        filestring = f.read()
    if not filestring.endswith('\n'):
        logger.error(f'{filename} must end with a newline character')
        raise FileFormatError(f'{filename} does not end with a newline')

def set_active_system_user_defined_gene(user_gene: dict, system_name:str) -> dict:
    user_gene['active_system'] = system_name
    user_gene['fasta']['active_system'] = system_name
    for model in user_gene['gene_model'].values():
        model['active_system'] = system_name
    return user_gene

def load_user_defined_gene(registry_file: Union[str, bytes, os.PathLike], system_name: str) -> UserDefinedGene:
    """
    Load a UserDefinedGene object from the registry
    
    This function is required because the active_system needs to be set for the schema prior to validation
    or else it will attempt to load files on alternate mount points that aren't reachable from the system in
    use.
    """
    try:
        with open(registry_file, 'r') as f:
            user_gene = json.load(f)
        logger.info(f'loaded {registry_file}')
        user_gene = set_active_system_user_defined_gene(user_gene=user_gene, system_name=system_name)
        # user_gene['active_system'] = system_name
        # user_gene['fasta']['active_system'] = system_name
        # for model in user_gene['gene_model'].values():
        #     model['active_system'] = system_name
        gene = UserDefinedGene.parse_obj(user_gene)
        # gene.propagate_active_system(system_name)
        logger.info(f'parsed {gene.id} gene model from {Path(registry_file).resolve()}')
    except Exception as e:
        logger.exception(f'failed to load gene model from {Path(registry_file).resolve()}\n{e}')
        raise
    return gene

def add_new_basepath(attribute, split_keyword: str, basepath: Union[str, bytes, os.PathLike], 
                     system_name: str, verify: bool=True) -> None:
    """
    Helper function for adding a new mount point to an existing registry from within a Genome or
    UserDefinedGene class object. This is used by the classes when adding new genomes or genes to a 
    registry that contains multiple mountpoints and is distinct from add_new_system_path() that is used
    when registering a new mountpoint.
    Splits a directory path on a keyword, appends a new basepath to create a new absolute
    path, and adds the new entry to the path dictionary using system_name as the key.
    """
    _, relpath = Path(next(iter(attribute.path.values()))).as_posix().split(f'{split_keyword}/')
    new_path = Path(os.path.join(basepath, split_keyword, relpath))
    try:
        if verify:
            if not new_path.exists():
                logger.error(f'file not found: {new_path}')
                raise FileNotFoundError(new_path)
            else:
                logger.debug(f'new filepath {new_path} is reachable')
        else:
            logger.debug(
                f'add_new_basepath without validation: ' +
                f'attribute: {attribute}; system_name: {system_name}l path: {basepath}')
    except Exception as e:
        logger.exception(e)
        raise

    if system_name in attribute.path:
        logger.info(
            f'system_name {system_name} is already a key in path -- skipping mountpoint update for {attribute}')
    else:
        attribute.path[system_name] = new_path
        logging.debug(f'{attribute} paths updated with "{system_name}": {new_path}')

def add_new_system_path(paths: dict, new_basepath: Union[str, bytes, os.PathLike], system_name: str,
                        split_keyword: str) -> dict[str, str]:
    _, relpath = Path(next(iter(paths.values()))).as_posix().split(f'{split_keyword}/')
    new_path = Path(os.path.join(new_basepath, split_keyword, relpath))
    paths[system_name] = new_path
    return paths

def add_new_genome_mountpoint(genome_dict: dict, mountpoint: Union[str, bytes, os.PathLike], system_name: str,
        verify: bool=True) -> dict:
    """Add a new mountpoint to the path dictionary of all class files with that attribute"""
    # logger.info(f'add_new_mountpoint called: mountpoint={mountpoint}; system_name={system_name}')
    for attribute in GENOME_CLASS_BASE_ATTR:
        genome_dict['base'][attribute]['path'] = add_new_system_path(
            paths=genome_dict['base'][attribute]['path'],
            new_basepath=mountpoint,
            system_name=system_name,
            split_keyword=GENOMES_RELATIVE_PATH
        )
        # add_new_basepath(model, attribute, GENOMES_RELATIVE_PATH,
        #                     basepath=mountpoint, system_name=system_name,
        #                     verify=verify)
    for attribute in GENOME_CLASS_MAIN_ATTR:
        genome_dict[attribute]['path'] = add_new_system_path(
            paths=genome_dict[attribute]['path'],
            new_basepath=mountpoint,
            system_name=system_name,
            split_keyword=GENOMES_RELATIVE_PATH
        )
        # add_new_basepath(model, attribute, GENOMES_RELATIVE_PATH,
        #                     basepath=mountpoint, system_name=system_name,
        #                     verify=verify)
    # logger.debug(f'add_new_mountpoint succesful')
    return genome_dict

def remove_genome_mountpoint(genome_dict: dict, system_name: str) -> dict:
    for attribute in GENOME_CLASS_BASE_ATTR:
        genome_dict['base'][attribute]['path'].pop(system_name)
    for attribute in GENOME_CLASS_MAIN_ATTR:
        genome_dict[attribute]['path'].pop(system_name)
    return genome_dict

def add_new_usergene_mountpoint(gene_dict: dict, mountpoint: Union[str, bytes, os.PathLike],
                                system_name: str) -> dict:
    """
    TODO
    """
    gene_dict['fasta']['path'] = add_new_system_path(
        paths=gene_dict['fasta']['path'],
        new_basepath=mountpoint,
        system_name=system_name,
        split_keyword=USER_GENES_RELATIVE_PATH)
    for version, in gene_dict['gene_model']:
        gene_dict['gene_model'][version]['path'] = add_new_system_path(
            paths=gene_dict['gene_model'][version]['path'],
            new_basepath=mountpoint,
            system_name=system_name,
            split_keyword=USER_GENES_RELATIVE_PATH
        )
    return gene_dict

def remove_usergene_mountpoint(gene_dict: dict, system_name: str) -> dict:
    """
    TODO
    """
    gene_dict['fasta']['path'].pop(system_name)
    for gene_model in gene_dict['gene_model'].values():
        gene_model['path'].pop(system_name)
    return gene_dict

def update_config_mountpoint(registry_path: Union[str, bytes, os.PathLike], system_name: str,
                             mode: str, remove_system_name: str=None) -> tuple[dict, dict]:
    """
    Iterates through all config files in a registry adding a new mount point entry to all
    paths using add_new_basepath. Returns a tuple of dictionaries. The 0-index holds genome
    configs and the 1-index holds user-defined-gene configs. The dictionary keys are paths to
    the config files and the values are the updated GenomeCollection objects that will replace
    them.
    """
    VALID_MODES = ['add', 'remove']
    if mode not in VALID_MODES:
        raise ValueError(f'mode must be one of {VALID_MODES}')
    # if mode == 'add' and mountpoint is None:
    #     raise ValueError(f"mountpoint is a required argument when mode='add'")
    if mode == 'remove' and remove_system_name is None:
        raise ValueError(f"remove_system_name is a required argument when mode='remove'")

    # update genome file paths
    genomes = {}
    user_defined_genes = {}
    for conf in os.listdir(Path(registry_path, GENOMES_CONFIG_DIR_RELATIVE_PATH).resolve()):
        conf_filepath = os.path.abspath(os.path.join(registry_path, GENOMES_CONFIG_DIR_RELATIVE_PATH, conf))
        # version_conf = load_genome(conf_filepath, system_name)
        with open(conf_filepath, 'r') as f:
            genome_collection = json.load(f)
        for genome_name, genome_dict in genome_collection['genomes'].items():
            logger.info(f'Updating {genome_name} (mode={mode})...')
            if mode == 'add':
                genome_collection['genomes'][genome_name] = add_new_genome_mountpoint(
                    genome_dict=genome_dict,
                    mountpoint=registry_path,
                    system_name=system_name
                )
                # genome.add_new_mountpoint(mountpoint=registry_path, system_name=system_name, verify=verify)
            elif mode == 'remove':
                genome_collection['genomes'][genome_name] =remove_genome_mountpoint(
                    genome_dict=genome_dict, system_name=remove_system_name)
        genomes[conf_filepath] = genome_collection

    # update user-defined-gene file paths
    for conf in os.listdir(Path(registry_path, USER_GENES_CONFIG_DIR_RELATIVE_PATH).resolve()):
        conf_filepath = os.path.abspath(os.path.join(registry_path, USER_GENES_CONFIG_DIR_RELATIVE_PATH, conf))
        # gene = load_user_defined_gene(conf_filepath, system_name)
        with open(conf_filepath, 'r') as f:
            gene_dict = json.load(f)
        logger.info(f"Updating {gene_dict['id']} (mode={mode})...")
        if mode == 'add':
            gene_dict = add_new_usergene_mountpoint(
                gene_dict=gene_dict,
                mountpoint=registry_path,
                system_name=system_name
            )
            # gene.add_new_mountpoint(mountpoint=registry_path, system_name=system_name)
        elif mode == 'remove':
            gene_dict = remove_usergene_mountpoint(gene_dict=gene_dict, system_name=remove_system_name)
            # gene.remove_mountpoint(remove_system_name=remove_system_name)
        user_defined_genes[conf_filepath] = gene_dict

    return genomes, user_defined_genes

def fetch_ensembl(url: str, destination_dir: Union[str, bytes, os.PathLike]) -> None:
    """
    Wrapper function to download a file from the Ensembl FTP server.
    Wraps urllib.request.urlopen and shutil.copyfileobj
    """
    target_filename = url.split('/')[-1]
    destination_filename = Path(destination_dir, target_filename)
    try:
        with urllib.request.urlopen(url) as response:
            logger.info(f"downloading {target_filename}: {humansize(response.getheader('content-length'))}")
            with open(destination_filename, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
                logger.info(f"download complete - {target_filename} stored as {out_file.name}")
    except HTTPError as e:
        logger.exception(f'{e.code}: {e.reason}; {url}\n{e}')
        if e.code == 404:
            logger.info(
                f'Tip: Check that assembly name matches the one available in the requested release. '
                f'See https://www.ensembl.org/info/website/archives/assembly.html to determine which assembly is matched to a release.')
        raise        


## Builder functions ##
def build_new_genome(registry_path: Union[str, bytes, os.PathLike], system_name: str,
                     genome_metadata: dict, input_dir: Union[str, bytes, os.PathLike]) -> Genome:
    """
    Create a new Genome object using a dictionary of metdata matching GenomeMetadata schema and
    and input_dir containing all required files to build a new Genome object: genome_fasta,
    gtf, transcriptome_fasta, star_index, refflat, and rrna.
    """
    # build registry directory name for this release and species
    genome_dir = Path(registry_path, GENOMES_RELATIVE_PATH, f"release-{genome_metadata['release']}",
        f"{genome_metadata['species']}")
    if genome_dir.exists():
        logger.error(f'a genome has already been stored in {genome_dir.resolve()}')
        raise FileExistsError(f'a genome has already been stored in {genome_dir.resolve()}')
    
    # get mount points
    mount_config = load_mount_config(registry_path)

    # get file paths from input directory
    try:
        genome_files = glob_genome_files(input_dir)
    except Exception as e:
        logger.exception(f"build new genome failed due to missing genome files\n{e}")
        raise

    try:
        write_paths = write_genome_files(
            genome_dir=genome_dir,
            genome_fasta=Path(genome_files['genome_fasta']),
            gtf=Path(genome_files['gtf']),
            transcriptome_fasta=Path(genome_files['transcriptome_fasta']),
            refflat=Path(genome_files['refflat']),
            rrna_interval_list=Path(genome_files['rrna_interval_list']),
            star_index=Path(genome_files['star_index'])
        )
    except FileExistsError as e:
        logger.exception(f'aborting due to attempted overwrite of existing registry files\n{e}')
        raise

    try:
        # the paths returned in write_paths are already absolute paths
        genome_fasta = GenomeFile(type='fasta', source='genome',
            default_system=mount_config.default_system_name,
            active_system=system_name,
            path={system_name: write_paths['genome_fasta']})
        gtf = GenomeFile(type='gtf',
            default_system=mount_config.default_system_name,
            active_system=system_name,
            path={system_name: write_paths['gtf']})
        transcriptome_fasta = GenomeFile(type='fasta', source='transcriptome',
            default_system=mount_config.default_system_name,
            active_system=system_name,
            path={system_name: write_paths['transcriptome_fasta']})
        star_index = GenomePath(type='star_index',
            default_system=mount_config.default_system_name,
            active_system=system_name,
            path={system_name: write_paths['star_index']})
        refflat = GenomeFile(type='refflat',
            default_system=mount_config.default_system_name,
            active_system=system_name,
            path={system_name: write_paths['refflat']})
        rrna_interval_list = GenomeFile(type='rrna_interval_list',
            default_system=mount_config.default_system_name,
            active_system=system_name,
            path={system_name: write_paths['rrna_interval_list']})
        base_genome = BaseGenome(
            metadata=genome_metadata,
            gtf=gtf,
            genome_fasta=genome_fasta
            )
        new_genome = Genome(
            id=base_genome.metadata.id,
            default_system=mount_config.default_system_name,
            base=base_genome,
            transcriptome_fasta=transcriptome_fasta,
            star_index=star_index,
            refflat=refflat,
            rrna_interval_list=rrna_interval_list,
            active_system=system_name)

        # add paths for any previously registered mountpoints
        for mount_sysname, mount_path in mount_config.mounts.items():
            new_genome.add_new_mountpoint(mountpoint=mount_path, system_name=mount_sysname, verify=False)

        conf_filename = f"{str(genome_metadata['release']).zfill(3)}.json"
        registry_file = Path(registry_path, GENOMES_CONFIG_DIR_RELATIVE_PATH, conf_filename)
        if registry_file.exists():
            genome_model = load_genome(registry_file, system_name)
            if genome_metadata['id'] in genome_model.genomes.keys():
                logger.error(f"aborting - genome with id {genome_metadata['id']} was already registered")
                raise DuplicateGenomeError(f"there is already a genome with id: {genome_metadata['id']}")
            else:
                genomes_list = genome_model.genomes
        else:
            genomes_list = dict()

        genomes_list[new_genome.id] = new_genome
        with registry_file.open('w') as f:
            f.write(GenomeCollection(genomes=genomes_list).json())
            logger.info(f'{getuser()} added genome {new_genome.id} to registry')

        logger.info(f'successfully built new genome {new_genome.id}')
        return new_genome
    except Exception as e:
        logger.exception(f'Exception raised in build_new_genome(); starting ERROR RECOVERY\n{e}')
        # try to delete any files and directories that were successfully created
        try:
            shutil.rmtree(genome_dir)
            logger.info(f'ERROR RECOVERY: {genome_dir} successfully removed')
        except FileNotFoundError:
            logger.info(f'ERROR RECOVERY: {genome_dir} was not created - no action required')
        # if this is the only species in this release, remove the release directory too
        if genome_dir.parent.is_dir():
            try:
                genome_dir.parent.rmdir()
                logger.info(f'ERROR RECOVERY: removed empty release directory {genome_dir.parent}')
            except Exception as e:
                logger.exception(f'failed to remove {genome_dir.parent}\n{e}')
                raise
        raise

def build_new_user_defined_gene(fasta: Union[str, bytes, os.PathLike], yaml_file: Union[str, bytes, os.PathLike],
        system_name: str, registry_path: Union[str, bytes, os.PathLike], **kwargs) -> UserDefinedGene:
    """
    Create a UserDefinedGene object from fasta file and YAML gene model for a single gene.
    If trying to add a collection of genes (e.g., ERCC) as a multi-fasta and associated YAML
    gene model definitions see build_new_user_defined_gene_collection().
    """
    try:
        model = gene_model_from_yaml(yaml_file)
        logger.info(model)
    except Exception as e:
        logger.exception(f'register-gene failed while loading gene model from YAML\n{e}')
        raise
    # always use absolute path in registry entries
    registry_path = Path(registry_path).resolve()
    target_dir = Path(registry_path, USER_GENES_RELATIVE_PATH, model.gene_id)
    if target_dir.exists():
        logger.error(
            f'{model.gene_id} has already been added to registry ({target_dir}), '
            f'try using update-gene if you are attempting to update the gene model for an already registered user-defined gene')
        raise FileExistsError(f'{target_dir} already exists')
    fasta_dest = Path(target_dir, f'{model.gene_id}.fa')
    yaml_fname = f'{model.gene_id}_v01.yaml'
    yaml_dest = Path(target_dir, yaml_fname)
    mount_config = load_mount_config(registry_path)
    try:
        target_dir.mkdir()
        shutil.copy(fasta, fasta_dest)
        shutil.copy(yaml_file, yaml_dest)
        yaml_file = GenomeFile(
            default_system=mount_config.default_system_name,
            active_system=system_name,
            path={system_name: yaml_dest}, type='yaml_gene_model')
        fasta_file = GenomeFile(
            default_system=mount_config.default_system_name,
            active_system=system_name, path={system_name: fasta_dest}, type='fasta')
        new_gene = UserDefinedGene(id=model.gene_id,
            default_system=mount_config.default_system_name,
            active_system=system_name, fasta=fasta_file, 
            gene_model={1: yaml_file})
        logger.info(model.gene_id)
        new_gene.propagate_active_system(system_name)

        # add paths for any previously registered mountpoints
        for mount_sysname, mount_path in mount_config.mounts.items():
            new_gene.add_new_mountpoint(mountpoint=mount_path, system_name=mount_sysname, verify=False)

        logger.info(f'built new user-defined-gene {model.gene_id}')

    except Exception as e:
        logger.exception(f'failed to build gene for fasta: {fasta}, YAML: {yaml_file}\n{e}')
        logger.info('starting error recovery')
        try:
            shutil.rmtree(target_dir)
            logger.info(f'ERROR RECOVERY: {target_dir} successfully removed')
        except FileNotFoundError:
            logger.info(f'ERROR RECOVERY: {target_dir} was not created - no action required')
            raise
        raise
    return new_gene


## Command line functions
def initialize(registry_path: Union[str, bytes, os.PathLike], system_name: str, group_name: Optional[str]=None, **kwargs) -> None:
    """
    Initialize a new genome registry
    
    Called via command line by `init`
    """

    # separate directories into lists of "open" and "restricted" to enable different permission settings
    # (even though current implementation gives the same permissions to both)
    open_subdirs = [USER_GENES_RELATIVE_PATH, USER_GENES_CONFIG_DIR_RELATIVE_PATH, LOG_DIR_RELATIVE_PATH]
    restricted_subdirs = [GENOMES_RELATIVE_PATH, GENOMES_CONFIG_DIR_RELATIVE_PATH]

    old_umask = os.umask(0o000)
    try:
        registry_path = Path(registry_path)
        if registry_path.is_file():
            raise FileExistsError
        elif registry_path.is_dir():
            if any(registry_path.iterdir()):
                raise FileExistsError(f'{registry_path.name} exists and is a non-empty directory')
            else:
                registry_path.chmod(0o775)
        else:
            registry_path.mkdir(0o775)
        for dirname in open_subdirs:
            new_path = Path(registry_path, dirname)
            for parent in reversed(new_path.parents):
                if not parent.exists():
                    parent.mkdir(0o775, exist_ok=True)
            new_path.mkdir(0o775)
        for dirname in restricted_subdirs:
            new_path = Path(registry_path, dirname)
            for parent in reversed(new_path.parents):
                if not parent.exists():
                    parent.mkdir(0o775, exist_ok=True)
            new_path.mkdir(0o775)

        # make log files and set permissions
        main_logfile = Path(registry_path, LOG_DIR_RELATIVE_PATH, 'genome-manager.log')
        main_logfile.touch()
        main_logfile.chmod(0o664)

        user_gene_logfile = Path(registry_path, LOG_DIR_RELATIVE_PATH, 'get-genes.log')
        user_gene_logfile.touch()
        user_gene_logfile.chmod(0o666)

        mount = {system_name: Path(registry_path)}
        mount_config = Mountpoints(default_system_name=system_name, mounts=mount)
        write_mount_config(registry_path=registry_path, mount_config=mount_config)
        Path(registry_path, MOUNTS_CONFIG_RELATIVE_PATH).chmod(0o664)
        logger = start_logger(registry_path=registry_path, command='_initialize')
        logger.info(f'Created new genome-registry: path={registry_path}, system_name={system_name}')

        if group_name:
            uid = os.getuid()
            gid = grp.getgrnam(group_name).gr_gid
            os.chown(registry_path, uid, gid)
            for root, dirs, files in os.walk(registry_path):
                for dir in dirs:
                    os.chown(Path(root, dir).resolve(), uid, gid)
                for file in files:
                    os.chown(Path(root, file).resolve(), uid, gid)
    # if a file or non-empty directory with this name exists DO NOT DELETE ANYTHING
    except FileExistsError:
        raise
    except:
        # undo everything
        if os.path.exists(registry_path) and os.path.isdir(registry_path):
            shutil.rmtree(registry_path)
        raise
    finally:
        os.umask(old_umask)

def register_genome(registry_path: Union[str, bytes, os.PathLike], system_name: str,
                    genome_metadata_file: Union[str, bytes, os.PathLike], input_dir: Union[str, bytes, os.PathLike],
                      **kwargs) -> None:
    """
    Add a new Genome to a GenomeCollection and write the json configuration file.
    Called via command line by `register-genome`.
    """
    registry_path = Path(registry_path).resolve()
    genome_metadata = parse_genome_metdata_file(genome_metadata_file)
    logger.info(f"{getuser()} called register-genome for {genome_metadata['id']}")

    try:
        build_new_genome(registry_path, system_name, genome_metadata, input_dir)
    except Exception as e:
        logger.exception(f'exception raised in register_genome(): {e}')
        raise

def delete_genome(genome):
    """Remove a genome from the registry"""
    raise NotImplementedError

def register_user_defined_gene(fasta: Union[str, bytes, os.PathLike], yaml_file: Union[str, bytes, os.PathLike],
        system_name: str, registry_path: Union[str, bytes, os.PathLike], **kwargs) -> None:
    """
    Add a new user-defined gene to the registry. Called by command line via `register-gene`
    See also update_user_defined_gene() to add a new gene model to an existing user-defined gene.
    """
    logger.info(f'{getuser()} called register-gene for registry: {Path(registry_path).resolve()}')
    mount_config = load_mount_config(registry_path)
    try:
        if not system_name in mount_config.mounts:
            raise InvalidSystemName
    except InvalidSystemName as e:
        logger.exception(
            f'failed to add new user-defined-gene - {system_name} is not a registered system-name for genome registry mounted at {registry_path}\n'
            f'check system-name spelling or use add-mountpoint to add a new mount point for this system first'
            f'\n{e}')
        raise
    except:
        raise

    validate_user_gene_file(fasta)
    validate_user_gene_file(yaml_file)

    new_gene = build_new_user_defined_gene(fasta=fasta, yaml_file=yaml_file, system_name=system_name,
        registry_path=registry_path)
    registry_file = Path(registry_path, USER_GENES_CONFIG_DIR_RELATIVE_PATH, new_gene.id + '.json')
    with registry_file.open('w') as f:
        f.write(new_gene.json())
    registry_file.chmod(0o775)
    logger.info(f'{getuser()} added user-defined gene {new_gene.id} version 1 to registry')

def update_user_defined_gene(registry_path: Union[str, bytes, os.PathLike],
        yaml_file: Union[str, bytes, os.PathLike], system_name: str, **kwargs) -> None:
    """
    Adds a new gene model to an existing user-defined gene that was added via `register-gene`
    Function is called by command line argument `update-gene`
    """
    logger.info(f'{getuser()} called update-gene for system_name: {system_name}; yaml_file: {yaml_file}')
    validate_user_gene_file(yaml_file)

    # load yaml_file and parse with YamlGeneModel
    try:
        model = gene_model_from_yaml(yaml_file)
    except Exception as e:
        logger.exception(f'update-gene failed while loading gene model from YAML\n{e}')
        raise
    # find the .json file
    # if it doesn't exist, then this isn't an update operation so fail with error message
    registry_path = Path(registry_path).resolve()
    registry_file = Path(registry_path, USER_GENES_CONFIG_DIR_RELATIVE_PATH, model.gene_id + '.json')
    if not registry_file.exists():
        logger.error(
            f'update-gene failed: no configuration file for {model.gene_id} found at {registry_file}')
        try:
            raise FileNotFoundError(f'{registry_file.resolve()} not found; try register-gene if adding gene for the first time')
        except Exception as e:
            logger.exception(e)
    gene = load_user_defined_gene(registry_file, system_name)

    # copy YAML to registry and update JSON config
    new_version_num, yaml_dest = gene.add_version(yaml_file, system_name)
    # UserDefinedGene.add_version() has error handling implemented so if it returns,
    # addition was successful and it's safe to update the config file.
    # However, an error writing the new config file would put the registry in a corrupted state.
    # To be safe, hold a copy of the original that can be restored and delete the added YAML file
    # on failure at this step.
    try:
        # hold onto a backup copy of original config in case it needs to be restored on failure
        with registry_file.open('r') as f:
            original_registry_file = f.read()
        with registry_file.open('w') as f:
            f.write(gene.json())
            logger.info(f'{getuser()} successfully updated {gene.id} to version {new_version_num}')
    except Exception as e:
        logger.exception(f'update-gene for {model.gene_id} encountered error updating the registry config file\n{e}')
        logger.info(f'ERROR RECOVERY: restoring previous config file version')
        logger.info(f'{model.gene_id} config JSON provided here as failsafe')
        logger.info(f'{original_registry_file}')
        with registry_file.open('w') as f:
            f.write(original_registry_file)
        if Path(yaml_dest).exists():
            logger.info(f'ERROR RECOVERY: deleting {yaml_dest} from registry')
            Path(yaml_dest).unlink()
        else:
            logger.info(f'ERROR RECOVERY: {yaml_file} was not written to registy - no further action required')

def get_user_defined_genes(registry_path: Union[str, bytes, os.PathLike], gene_ids: Union[str, list[str]],
        system_name: str, outdir: Union[str, bytes, os.PathLike] = Path('.'), 
        version_delim: str='.', **kwargs) -> None:
    """
    Retrieve user-defined gene(s) and write files to disk. If multiple genes, concatenate files
    before returning.
    Called from command line with command `get-genes`
    """
    if isinstance(gene_ids, str):
        gene_ids = [gene_ids]
    if not isinstance(gene_ids, list):
        raise TypeError('gene_ids must be of type str or list of strings')
    if not all([isinstance(gene, str) for gene in gene_ids]):
        raise TypeError('gene_ids must be a string or list of strings')

    fasta_files = []
    yaml_files = []
    collected_ids = []

    for gid in gene_ids:
        if version_delim in gid:
            gene_id, version = gid.strip().split(version_delim)
        else:
            gene_id = gid.strip()
            version = -1

        registry_file = Path(registry_path, USER_GENES_CONFIG_DIR_RELATIVE_PATH, gene_id + '.json')
        gene = load_user_defined_gene(registry_file, system_name)
        fasta_files.append(Path(gene.fasta.path[system_name]))
        yaml_files.append(gene.get_version(version, system_name))
        collected_ids.append(gene.id)
    if len(collected_ids) > 3:
        basename = 'custom'
    else:
        basename = '.'.join(collected_ids)

    with open(Path(outdir, basename + '.fa'), 'wb') as outfile:
        for f in fasta_files:
            with open(f, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)

    yaml_concat = fileinput.input(yaml_files)
    yaml_out = Path(outdir, basename + '.yaml')
    with open(yaml_out, 'w') as outfile:
        for line in yaml_concat:
            outfile.write(line)

    gene_model = generate_gtf_entry.read_genes_from_yaml(yaml_out)
    gtf = generate_gtf_entry.generate_gtf(gene_model)
    with open(Path(outdir, basename + '.gtf'), 'w') as outfile:
        outfile.write(gtf)

def add_mountpoint(registry_path: Union[str, bytes, os.PathLike], system_name: str, **kwargs) -> None:
    """
    Implements the add-mountpoint command which updates all config files with new absolute paths
    for a different mount point, keyed by system_name. No genome or gene data files are added
    or modified via this command. The intended usage is to update absolute paths to access 
    genome registry files stored on a single NAS that might be mounted to different paths on
    different systems. All new filepaths are verified.

    Called from command line by command `add-mountpoint`
    """
    logger.info(f'add-mountpoint: registry-path={registry_path}, system-name={system_name}')
    try:
        if not Path(registry_path).is_dir():
            raise FileNotFoundError(
                f'{registry_path} is not a reachable path: add-mountpoint must be performed from the system to be addded')
    except Exception as e:
        logger.exception(e)
        raise

    try:
        mount_config = load_mount_config(registry_path)
        for sysname, mount_path in mount_config.mounts.items():
            if str(system_name) == str(sysname):
                raise ValueError(
                    f'system_name {system_name} has already been used to register mountpoint {mount_path}')
            if Path(registry_path) == Path(mount_path):
                raise ValueError(
                    f'registry-path {registry_path} has already been added as system-name {sysname}')
    except Exception as e:
        logger.exception(e)
        raise

    try:
        genome_conf_directory = Path(registry_path, GENOMES_CONFIG_DIR_RELATIVE_PATH)
        genome_recovery_dir = shutil.copytree(genome_conf_directory,
                                    Path(registry_path, TEMP_DIR_RELATIVE_PATH, 'genome_recovery'),
                                    dirs_exist_ok=True)
        user_conf_directory = Path(registry_path, USER_GENES_CONFIG_DIR_RELATIVE_PATH)
        user_recovery_dir = shutil.copytree(user_conf_directory,
                                    Path(registry_path, TEMP_DIR_RELATIVE_PATH, 'user_gene_recovery'),
                                    dirs_exist_ok=True)
    except Exception as e:
        # do not leave behind any backup directories if try block failed
        logger.exception(e)
        if 'genome_recovery_dir' in locals() and genome_recovery_dir.exists():
            logger.info(f'Removing backup: {genome_recovery_dir}')
            shutil.rmtree(genome_recovery_dir)
        if 'user_recovery_dir' in locals() and user_recovery_dir.exists():
            logger.info(f'Removing backup: {user_recovery_dir}')
            shutil.rmtree(user_recovery_dir)
        raise

    try:
        new_genomes, new_genes = update_config_mountpoint(
            registry_path=registry_path, system_name=system_name, mode='add')
        for config, genome_dict in new_genomes.items():
            logger.info(f'updating mountpoint for {config}')
            genome_dict = set_active_system_genome(genome_dict, system_name)
            genome = GenomeCollection(**genome_dict)
            with open(config, 'w') as f:
                f.write(genome.json())
        for config, gene_dict in new_genes.items():
            logger.info(f'updating mountpoint for {config}')
            gene_dict = set_active_system_user_defined_gene(gene_dict, system_name)
            gene = UserDefinedGene(**gene_dict)
            with open(config, 'w') as f:
                f.write(gene.json())
        mount_config.mounts[system_name] = Path(registry_path)
        write_mount_config(registry_path=registry_path, mount_config=mount_config)
    except Exception as e:
        logger.exception(e)
        logger.info('mount point addition failed -- restoring old config')
        shutil.copytree(genome_recovery_dir, genome_conf_directory, dirs_exist_ok=True)
        shutil.copytree(user_recovery_dir, user_conf_directory, dirs_exist_ok=True)
        raise
    finally:
        logger.info('removing temporary recovery directories (if any)')
        if genome_recovery_dir.exists():
            shutil.rmtree(genome_recovery_dir)
            logger.info('removed genome recovery directory')
        if user_recovery_dir.exists():
            shutil.rmtree(user_recovery_dir)
            logger.info('removing user-defined-genes recovery directory')

def remove_mountpoint(registry_path: Union[str, bytes, os.PathLike], remove_system_name: str,
                      **kwargs) -> None:
    """
    Remove a previously registered mountpoint indexed by remove_system_name. The system_name argument
    defines the system that the command is being issued from while remove_system_name is the one to remove.

    Called via command line by `remove-mountpoint`
    """
    active_system_name = find_active_system(registry_path)
    logger.info(f'remove-mountpoint: registry-path={registry_path}, system-name={remove_system_name}, call from system-name={active_system_name}')
    try:
        mount_config = load_mount_config(registry_path)
        if str(remove_system_name) == str(active_system_name):
            raise ValueError(
                f"Cannot remove system_name {remove_system_name} because it is the active one."
                f"The remove-mountpoint command must be issued from a registered system-name that is not the one being removed.")
        if mount_config.default_system_name == remove_system_name:
            raise ValueError(f'Cannot remove default_system: {remove_system_name} is the default_system for the registry')
    except Exception as e:
        logger.exception(e)
        raise

    confirmation = input(
        f"\nThis action will permamently remove the following mount point from the registry:\n"
        f"\tsystem-name: {remove_system_name}\n"
        f"\tmount path: {mount_config.mounts[remove_system_name]}\n\n"
        f"Type CONFIRM to confirm this action: ")
    if confirmation != 'CONFIRM':
        logger.info(f"remove-mountpoint canceled by user {getuser()}")
        return
    else:
        logger.info(f"remove-mountpoint confirmed by user {getuser()}")

    try:
        genome_conf_directory = Path(registry_path, GENOMES_CONFIG_DIR_RELATIVE_PATH)
        genome_recovery_dir = shutil.copytree(genome_conf_directory,
                                    Path(registry_path, TEMP_DIR_RELATIVE_PATH, 'genome_recovery'),
                                    dirs_exist_ok=True)
        user_conf_directory = Path(registry_path, USER_GENES_CONFIG_DIR_RELATIVE_PATH)
        user_recovery_dir = shutil.copytree(user_conf_directory,
                                    Path(registry_path, TEMP_DIR_RELATIVE_PATH, 'user_gene_recovery'),
                                    dirs_exist_ok=True)
    except Exception as e:
        # do not leave behind any backup directories if try block failed
        logger.exception(e)
        if 'genome_recovery_dir' in locals() and genome_recovery_dir.exists():
            logger.info(f'Removing backup: {genome_recovery_dir}')
            shutil.rmtree(genome_recovery_dir)
        if 'user_recovery_dir' in locals() and user_recovery_dir.exists():
            logger.info(f'Removing backup: {user_recovery_dir}')
            shutil.rmtree(user_recovery_dir)
        raise

    try:
        new_genomes, new_genes = update_config_mountpoint(
            registry_path=registry_path, system_name=active_system_name,
            mode='remove', remove_system_name=remove_system_name)
        for config, genome in new_genomes.items():
            logger.info(f'updating mountpoint for {config}')
            genome = set_active_system_genome(genome, active_system_name)
            genome_collection = GenomeCollection(**genome)
            with open(config, 'w') as f:
                f.write(genome_collection.json())
        for config, gene in new_genes.items():
            logger.info(f'updating mountpoint for {config}')
            gene = set_active_system_user_defined_gene(gene, active_system_name)
            user_defined_gene = UserDefinedGene(**gene)
            with open(config, 'w') as f:
                f.write(user_defined_gene.json())
        mount_config.mounts.pop(remove_system_name, None)
        write_mount_config(registry_path=registry_path, mount_config=mount_config)
    except Exception as e:
        logger.exception(e)
        logger.info('mount point removal failed -- restoring old config')
        shutil.copytree(genome_recovery_dir, genome_conf_directory, dirs_exist_ok=True)
        shutil.copytree(user_recovery_dir, user_conf_directory, dirs_exist_ok=True)
        raise
    finally:
        logger.info('removing temporary recovery directories (if any)')
        if genome_recovery_dir.exists():
            shutil.rmtree(genome_recovery_dir)
            logger.info('removed genome recovery directory')
        if user_recovery_dir.exists():
            shutil.rmtree(user_recovery_dir)
            logger.info('removing user-defined-genes recovery directory')

def list_genomes(registry_path: Union[str, bytes, os.PathLike], system_name: str=None,
                 species: str=None, **kwargs) -> str:
    """
    List all the currently registered genomes listing the ID to use when calling it as well the 
    assembly and release. If optional argument 'species' is not provided, all species will be listed
    in alphabetical order. Otherwise, only genomes for the selected species will be displayed.
    
    Called via command line by `list-genomes`
    """
    genomes_combined = defaultdict(list)
    # use the mountpoint config file to find the system name if it is not provided
    if system_name is None:
        try:
            system_name = find_active_system(registry_path)
            # system_found = False
            # mount_config = load_mount_config(registry_path)
            # for sysname, regpath in mount_config.mounts.items():
            #     if Path(regpath) == Path(registry_path):
            #         system_name = sysname
            #         system_found = True
            #         break
            # if not system_found:
            #     # this should not be possible, but include for debugging purposes just in case
            #     raise ValueError(f'could not find system_name for registry_path: {registry_path}')
        except Exception as e:
            logger.exception(e)
            raise

    for filename in Path(registry_path, GENOMES_CONFIG_DIR_RELATIVE_PATH).glob('*.json'):
        config = load_genome(filename, system_name)
        for current_species, metadata in config.get_genome_info().items():
            genomes_combined[current_species].append(metadata)

    for val in genomes_combined.values():
        val.sort(key=lambda x: int(x.get('release')))
    sorted_genomes = dict(sorted(genomes_combined.items()))

    genomes_found = False
    for current_species, releases in sorted_genomes.items():
        if species is None or current_species == species:
            genomes_found = True
            print(f'{current_species}:')
            for release in releases:
                print(f'{release.get("id"): <14}{release.get("assembly") + ", release " + str(release.get("release")): <20}')
            print()

    if not genomes_found:
        species_str = f" for species '{species}'" if species is not None else ''
        print(f'No registered genomes found{species_str}.')

def list_user_defined_genes(registry_path: Union[str, bytes, os.PathLike], system_name: str=None,
                            **kwargs) -> str:
    """
    List all the currently registered user-defined gene IDs.
    
    Called via command line using `list-genes`
    """
    # mount_config = load_mount_config(registry_path)
    # use the mountpoint config file to find the system name if it is not provided
    if system_name is None:
        try:
            system_name = find_active_system(registry_path)
            # system_found = False
            # mount_config = load_mount_config(registry_path)
            # for sysname, regpath in mount_config.mounts.items():
            #     if Path(regpath) == Path(registry_path):
            #         system_name = sysname
            #         system_found = True
            #         break
            # if not system_found:
            #     # this should not be possible, but include for debugging purposes just in case
            #     raise ValueError(f'could not find system_name for registry_path: {registry_path}')
        except Exception as e:
            logger.exception(e)
            raise
    user_defined_genes = []
    for config_file in Path(registry_path, USER_GENES_CONFIG_DIR_RELATIVE_PATH).glob('*.json'):
        gene = load_user_defined_gene(config_file, system_name)
        user_defined_genes.append(gene.id)
    print('Available user-defined genes by ID:')
    print('\n'.join(sorted(user_defined_genes)))

def list_mountpoints(registry_path: Union[str, bytes, os.PathLike], **kwargs) -> str:
    """
    List all the currently registered mount points in a 2-column format of <system-name> <mount point>
    
    Called via command line using `list-mountpoints`
    """
    mount_config = load_mount_config(registry_path)
    print(mount_config)

def download_ensembl_genome(registry_path: Union[str, bytes, os.PathLike], species: str, release: int,
                            assembly_name: str=None, use_cwd: bool=False, **kwargs) -> GenomeMetadata:
    """
    Download a genome fasta and gtf file from Ensembl. These are the source files for genomes in the registry.
    By default, the files are downloaded to a temporary download directory within the genome registry, but this
    can be overridden by setting use_cwd=True, in which case the files are downloaded to the current working
    directory instead.

    Called from command line by command `download-genome`
    """
    ASSEMBLY_NAME_DICT = {
        'homo_sapiens': 'GRCh38',
        'mus_musculus': 'GRCm39',
        'rattus_norvegicus': 'mRatBN7.2',
        'macaca_fascicularis': 'Macaca_fascicularis_6.0',
        'macaca_mulatta': 'Mmul_10',
        'sus_scrofa': 'Sscrofa11.1',
        'cricetulus_griseus': 'CHOK1GS_HDv1',
        'monodon_monoceros': 'NGI_Narwhal_1'
    }

    try:
        assembly_name = ASSEMBLY_NAME_DICT[species.lower()] if assembly_name is None else assembly_name
    except KeyError as e:
        logger.info(f"Cannot find assembly name for species '{species}'."
                    f" Check species spelling or provide 'assembly' argument directly if species is spelled correctly.")
        logger.exception(e)
        raise
    gtf_target_filename = f'{species.capitalize()}.{assembly_name}.{release}.gtf.gz'
    gtf_url = f'https://ftp.ensembl.org/pub/release-{release}/gtf/{species.lower()}/{gtf_target_filename}'

    fasta_primary_target_filename = f'{species.capitalize()}.{assembly_name}.dna.primary_assembly.fa.gz'
    fasta_primary_url = f'https://ftp.ensembl.org/pub/release-{release}/fasta/{species.lower()}/dna/{fasta_primary_target_filename}'

    fasta_toplevel_target_filename = f'{species.capitalize()}.{assembly_name}.dna.toplevel.fa.gz'
    fasta_toplevel_url = f'https://ftp.ensembl.org/pub/release-{release}/fasta/{species.lower()}/dna/{fasta_toplevel_target_filename}'

    destination_dir = Path(os.getcwd()) if use_cwd else Path(
        registry_path, TEMP_DOWNLOAD_RELATIVE_PATH, f'release-{release}', species.lower())
    if not destination_dir.exists():
        destination_dir.mkdir(parents=True, exist_ok=True)

    # get gtf file
    fetch_ensembl(gtf_url, destination_dir)

    # get fasta file: try for primary_assembly first, but if it doesn't exist, then toplevel is the same thing as documented in Ensembl readme files
    try:
        fetch_ensembl(fasta_primary_url, destination_dir)
        assembly_type = 'primary_assembly'
    except HTTPError:
        fetch_ensembl(fasta_toplevel_url, destination_dir)
        assembly_type = 'toplevel'

    metadata = GenomeMetadata(
        id=f"{format_assembly_name(assembly_name)}:{release}",
        species=species,
        species_short=abbreviate_species(species),
        release=release,
        assembly=assembly_name,
        assembly_type=assembly_type,
        sequence_type='dna')
    
    with open(Path(destination_dir, 'metadata.json'), 'w') as f:
        f.write(metadata.json())

    return metadata

def clean(registry_path: Union[str, bytes, os.PathLike], **kwargs) -> None:
    """
    Function to clean the registry temp directories, if any exist

    Called from command line by command `clean`
    """
    target = Path(registry_path, TEMP_DIR_RELATIVE_PATH).resolve()
    if target.exists():
        total_size = 0
        num_files = 0
        for dirpath, dirnames, filenames in os.walk(target):
            for f in filenames:
                fp = os.path.join(dirpath, f)
                # skip if it is symbolic link
                if not os.path.islink(fp):
                    total_size += os.path.getsize(fp)
                    num_files += 1
        logging.info(f"found {num_files} temporary files totaling {humansize(total_size)}")
        shutil.rmtree(target)
        logging.info(f"removed all temporary directories and files recovering {humansize(total_size)} of space")
    else:
        logging.info("no temporary directories or files found")


## Logger/Parser/Main ##
def start_logger(registry_path: Union[str, bytes, os.PathLike], **kwargs) -> None:
    if kwargs['command'] == 'init':
        # no logging for init
        return
    elif kwargs['command'] == 'get-genes':
        logfile = Path(registry_path, LOG_DIR_RELATIVE_PATH, 'get-genes.log')
        logfile.touch()
        logfile.chmod(0o777)
    else:
        logfile = Path(registry_path, LOG_DIR_RELATIVE_PATH, 'genome-manager.log')

    try:
        logging.basicConfig(
            filename=logfile,
            encoding='utf-8',
            level=logging.DEBUG,
            format='%(asctime)s [%(levelname)s]: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        logger = logging.getLogger()
    except FileNotFoundError as e:
        raise FileNotFoundError(f'{e}\ncheck that registry-path points to a valid genome-registry')

    # log to sys.stderr as well, unless it is a list-* command in which case we want the output to be clean
    if kwargs['command'] not in ['list-mountpoints', 'list-genomes', 'list-genes']:
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.INFO)
        logger.addHandler(stream_handler)
    logger.info(
        f'genome_manager.py v{__version__} called by {getuser()} using Python {python_version()}')
    logger.info(f'command line: {" ".join(quote(s) for s in sys.argv)}')
    return logger

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action='version',
        version='%(prog)s {version}'.format(version=__version__))
    sp = parser.add_subparsers(dest='command')

    init_parser = sp.add_parser('init', help='initialize a new genome registry')
    init_parser.set_defaults(func=initialize)
    init_parser.add_argument('--registry-path', required=True, help='path to the genome registry')
    init_parser.add_argument('--system-name', required=True,
        help = 'string identifying the system where paths will exist (e.g., HPC name, workstation name, etc.)')
    init_parser.add_argument('--group-name',
        help='optional permission group name for write access to non-restricted subdirectories (e.g., user_defined_genes)')

    register_genome_parser = sp.add_parser('register-genome', help='register a genome')
    register_genome_parser.set_defaults(func=register_genome)
    register_genome_parser.add_argument('--registry-path', required=True,
        help='path to the genome registry')
    register_genome_parser.add_argument('--system-name', required=True,
        help = 'string identifying the system where paths will exist (e.g., HPC name, workstation name, etc.)')
    register_genome_parser.add_argument('--genome-metadata-file', required=True,
        help = 'path to JSON file containing genome metadata')
    register_genome_parser.add_argument('--input-dir', required=True,
        help='path to input directory containing all required genome files')

    register_gene_parser = sp.add_parser('register-gene', help='register a user-defined gene')
    register_gene_parser.set_defaults(func=register_user_defined_gene)
    register_gene_parser.add_argument('--registry-path', required=True,
        help='path to the genome registry')
    register_gene_parser.add_argument('--yaml-file', required=True,
        help = 'path to YAML file containing gene model for custom fasta')
    register_gene_parser.add_argument('--fasta', required=True,
        help = 'path to fasta file containing sequence of a genome modification')
    register_gene_parser.add_argument('--system-name', required=True,
        help = 'string identifying the system where paths will exist (e.g., HPC name, workstation name, etc.)')

    update_gene_parser = sp.add_parser('update-gene', help='update an existing user-defined gene with a new YAML gene model')
    update_gene_parser.set_defaults(func=update_user_defined_gene)
    update_gene_parser.add_argument('--registry-path', required=True,
        help='path to the genome registry')
    update_gene_parser.add_argument('--yaml-file', required=True,
        help = 'path to YAML file containing gene model for custom fasta')
    update_gene_parser.add_argument('--system-name', required=True,
        help = 'string identifying the system where paths will exist (e.g., HPC name, workstation name, etc.)')

    get_genes_parser = sp.add_parser('get-genes',
        help='retrieve fasta and YAML genes models for selected user-defined genes')
    get_genes_parser.set_defaults(func=get_user_defined_genes)
    get_genes_parser.add_argument('--registry-path', required=True,
        help='path to the genome registry')
    get_genes_parser.add_argument('--gene-ids', required=True, nargs='+',
        help = 'comma-separated list of gene IDs to retrieve from registry')
    get_genes_parser.add_argument('--system-name', required=True,
        help = 'string identifying the system where paths will exist (e.g., HPC name, workstation name, etc.)')
    get_genes_parser.add_argument('--outdir', required=False, default='./',
        help = 'target directory to write output to')
    get_genes_parser.add_argument('--version-delim', required=False, default='.',
        help = "delimiter when specifying a specific gene model version number after a gene-id (default = '.')")

    add_mountpoint_parser = sp.add_parser('add-mountpoint',
                                          help='add absolute paths for a new mount point when a single physical genome registry is shared across a network')
    add_mountpoint_parser.set_defaults(func=add_mountpoint)
    add_mountpoint_parser.add_argument('--registry-path', required=True,
        help='path to be added to the genome registry (i.e., the top-level genome registry directory)')
    add_mountpoint_parser.add_argument('--system-name', required=True,
        help = 'string identifying the system where new registry-path exists (e.g., HPC name, workstation name, etc.)')

    remove_mountpoint_parser = sp.add_parser('remove-mountpoint',
                                          help='remove a previously registered mountpoint by system-name')
    remove_mountpoint_parser.set_defaults(func=remove_mountpoint)
    remove_mountpoint_parser.add_argument('--registry-path', required=True,
        help='path to the genome registry')
    remove_mountpoint_parser.add_argument('--remove-system-name', required=True,
        help = 'the system-name string of the mount point to be removed')

    list_mountpoints_parser = sp.add_parser('list-mountpoints',
        help='list all registered mount points in a two column format of <system-name> <mount point>')
    list_mountpoints_parser.set_defaults(func=list_mountpoints)
    list_mountpoints_parser.add_argument('--registry-path', required=True,
        help='path to the genome registry')

    list_genomes_parser = sp.add_parser('list-genomes',
        help='list all registered genomes arranged by species')
    list_genomes_parser.set_defaults(func=list_genomes)
    list_genomes_parser.add_argument('--registry-path', required=True,
        help='path to the genome registry')
    list_genomes_parser.add_argument('--species', required=False,
        help='optional species to search for')

    list_genes_parser = sp.add_parser('list-genes',
        help='list all registered user-defined genes by ID')
    list_genes_parser.set_defaults(func=list_user_defined_genes)
    list_genes_parser.add_argument('--registry-path', required=True,
        help='path to the genome registry')

    download_genome_parser = sp.add_parser('download-genome', help='download the source genome files from Ensembl FTP server')
    download_genome_parser.set_defaults(func=download_ensembl_genome)
    download_genome_parser.add_argument('--registry-path', required=True,
        help='path to the genome registry')
    download_genome_parser.add_argument('--species', required=True,
        help='full species name as a string with underscore, e.g., mus_musculus')
    download_genome_parser.add_argument('--release', required=True,
        help='Ensembl release number as a string or integer')
    download_genome_parser.add_argument('--use_cwd', action='store_true',
        help='boolean (default False): when True, download files to current working directory rather than genome registry')
    download_genome_parser.add_argument('--assembly-name', required=False,
        help='optional string stating the assembly name name, e.g., GRCh38; not necessary for latest builds of common species (human, mouse, rat, cyno)')

    clean_parser = sp.add_parser('clean', help='delete any temporary files/directories to recovery space')
    clean_parser.set_defaults(func=clean)
    clean_parser.add_argument('--registry-path', required=True,
        help='path to the genome registry')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_arguments()
    logger = start_logger(**vars(args))
    args.func(**vars(args))
else:
    logger = logging.getLogger(__name__)
    logger.addHandler(logging.StreamHandler())
