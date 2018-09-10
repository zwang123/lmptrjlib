from classproperty import classproperty, classproperty_support
from collections import OrderedDict

@classproperty_support
class LMPtrj(object):
    """
    A class that parses lammps trajectories and stores the data in the
    trj attribute as a dictionary

    attributes:
    clear() : clear the trj, called by the class constructor
    parse(trjname) : construct the trj, called by the class constructor
    to_file() : write trj to a file
    sort() : sort the atoms section
    sort_id() : sort the atoms section based on id

    valid_sections: the keys to the trj.values()
    data_types: atoms section, key to type dictionary
    trj     : a dictionary, key = timestep, value = dictionary

    {timestep : {valsection  : value, ...,
                 dictsection : (args, [{ attr: value}, ...]), ...}, ...}

    if two sections have the same timestep, the latter one will overwrite 
    the former one
    """
    def __init__(self, trjname=None):
        self.clear()
        self.parse(trjname)

########################## public methods ##########################

    def clear(self):
        self.__trj = OrderedDict()

    def parse(self, trjname):
        if trjname is None:
            return self

        with open(trjname) as f:
            data = f.readlines()

        next_item_idx = 0
        timestep = None
        while len(data) > next_item_idx:
            line_idx = next_item_idx
            line = data[line_idx]

            if line[:5] != "ITEM:":
                raise RuntimeError("Section not starting with ITEM")

            section, args = self._parse_item(line[5:].strip())
            section_value, lines_parsed = self.__getattribute__( 
                    "_parse_" + section)\
                    (args, timestep, data, line_idx + 1)
            next_item_idx = line_idx + lines_parsed + 1

            if section == "timestep":
                timestep = section_value
                self.__trj[timestep] = OrderedDict()

            self.__set_section_value(timestep, section, section_value)

        return self

    def to_file(self, filename=None):
        data = ""
        for timestep, frame in self.trj.items():
            for section, value in frame.items():
                try:
                    args, sectionlist = value
                except TypeError:
                    args = []
                    sectionlist = [{None: value}]
                data += "ITEM: {:s}\n".format(' '.join(
                    [self.__valid_keys_rev[section]] + args))
                for line in sectionlist:
                    data += " ".join([str(x) for x in line.values()]) + "\n"
        if filename:
            with open(filename, "w") as f:
                f.write(data)
        return data

    def sort(self, key):
        for frame in self.trj.values():
            frame["atoms"][1].sort(key=key)
        return self

    def sort_id(self):
        return self.sort(lambda x : x["id"])

    def sort_mol_type_id(self):
        return self.sort(lambda x : [x["mol"], x["type"], x["id"]])
                
########################## public members ##########################

    @property
    def trj(self):
        return self.__trj

    @classproperty
    def valid_sections(cls):
        return cls.__valid_sections

    @classproperty
    def data_types(cls):
        return cls.__data_types

########################## private methods ##########################

    def __set_section_value(self, timestep, section, value):
        self.__trj[timestep][section] = value

########################## protected methods ##########################

    def _parse_atoms(self, args, timestep, data, start):
        natoms = self.trj[timestep]["natoms"]
        assert(len(data) >= start+natoms)
        return (args, self._parse_formatted_section(args, 
                [self.data_types[arg] for arg in args],
                data[start:start+natoms])), natoms
    
    @classmethod
    def _parse_item(cls, line):
        for key in cls.__valid_sections:
            if key == line[:len(key)]:
                return cls.__valid_keys[key], line[len(key):].split()
        raise ValueError(
                'Does not recognize the item "{}". Valid sections are '.\
                        format(line) + ', '.join(cls.__valid_sections) + '.')
        return None, None

    @classmethod
    def _parse_bounds(cls, args, timestep, data, start):
        for arg in args:
            assert(arg == "pp") 
        ndim = len(args)
        assert(len(data) >= start+ndim)
        return (args, cls._parse_formatted_section(["lbound", "hbound"], 
                [float, float],
                data[start:start+ndim])), ndim

    @staticmethod
    def _parse_natoms(args, timestep, data, start):
        """ Must not use timestep, which is uninitialized, in this function """
        if args:
            raise RuntimeError("Non empty args")
        return int(data[start]), 1

    _parse_timestep = _parse_natoms

    @staticmethod
    def _parse_formatted_section(args, types, data):
        rtn = []
        for line in data:
            elem = OrderedDict()
            for attr, type_, val in zip(args, types, line.split()):
                elem[attr] = type_(val)
            rtn.append(elem)
        return rtn

########################## private static members ##########################

    # could be overridden by derived class
    __valid_keys = OrderedDict([
            ("TIMESTEP" , "timestep"), 
            ("NUMBER OF ATOMS", "natoms"), 
            ("BOX BOUNDS" , "bounds"), 
            ("ATOMS" , "atoms"),
            ])
    __valid_keys_rev = OrderedDict([
        (v, k) for k, v in __valid_keys.items()
            ])
    __valid_sections = __valid_keys.keys()
    __data_types = OrderedDict([
            ("id" , int),
            ("mol" , int),
            ("type" , int),
            ("q" , float),
            ("x" , float),
            ("y" , float),
            ("z" , float),
            ("vx" , float),
            ("vy" , float),
            ("vz" , float),
            ("fx" , float),
            ("fy" , float),
            ("fz" , float),
            ])

if __name__ == "__main__":
    foo = LMPtrj()
    bar = LMPtrj("test/test.lammpstrj")
    foo.parse("test/test.lammpstrj")
    foo.to_file("result/rtn.lammpstrj")
    assert(foo.trj == bar.trj)
    foo.sort_id().to_file("result/sorted_id.lammpstrj")
    foo.sort_mol_type_id().to_file("result/sorted_moltypeid.lammpstrj")
    #print(foo.trj)
    #print(foo.valid_sections)
