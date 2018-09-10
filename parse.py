#from util import constclassproperty
from classproperty import classproperty, classproperty_support

@classproperty_support
class LMPtrj(object):
    """
    A class that parses lammps trajectories and stores the data in the
    trj attribute as a dictionary

    attributes:

    trj     : a dictionary, key = timestep, value = dictionary

    {timestep : {valsection  : value, ...,
                 dictsection : [{ attr: value}, ...], ...}, ...}

    valid_sections: the keys to the trj.values()
    parse(trjname) : construct the trj, called by the class constructor
    """
    def __init__(self, trjname=None):
        self.clear()
        self.parse(trjname)

    #@staticmethod
    #def __foo(self):
    #    pass

    def clear(self):
        self.__trj = dict()


    def parse(self, trjname):
        #print(self.__getattribute__("__foo"))
        #self._parse_item("dfasd")
        #self.__trj = dict()

        if trjname is None:
            return self

        with open(trjname) as f:
            data = f.readlines()

        #for line_idx, line in enumerate(data):
        next_item_idx = 0
        timestep = None
        while len(data) > next_item_idx:

            line_idx = next_item_idx
            line = data[line_idx]

            #if line[:5] == "ITEM:":
            if line[:5] != "ITEM:":
                raise RuntimeError("Section not starting with ITEM")

            section, args = self._parse_item(line[5:].strip())

            #print(self.__getattribute__("parse"))
            #print(self.__getattribute__("data_types"))
            #print(self.__getattribute__("__set_section_value"))
            #print(self.__getattribute__("__trj"))
            #print(vars(type(self)))
            #exit()
            #section_value, lines_parsed = self.__dict__["_parse_" + section]\
            #section_value, lines_parsed = getattr(self, 
            section_value, lines_parsed = self.__getattribute__( 
                    "_parse_" + section)\
                    (args, timestep, data, line_idx + 1)
            next_item_idx = line_idx + lines_parsed + 1


            if section == "timestep":
                #if args:
                #    raise RuntimeError("non empty args for TIMESTEP")
                #timestep = int(data[line_idx + 1])
                #section_value = timestep
                timestep = section_value
                self.__trj[timestep] = dict()
            #else:
            self.__set_section_value(timestep, section, section_value)

            #if len(data) <= next_item_idx:
            #    break
            #elif data[next_item_idx][:5] != "ITEM:":
            #    raise RuntimeError("Unknown line after section {}".\
            #            format(section))
                
            # TODO
            #make sure
            #print(section, args)
        return self

    @property
    def trj(self):
        return self.__trj

    #@property
    #@constclassproperty
    @classproperty
    def valid_sections(cls):
        return cls.__valid_sections

    @classproperty
    def data_types(cls):
        return cls.__data_types

    def __set_section_value(self, timestep, section, value):
        self.__trj[timestep][section] = value
    
    @classmethod
    def _parse_item(cls, line):
        for key in cls.__valid_sections:
            if key == line[:len(key)]:
                return cls.__valid_keys[key], line[len(key):].split()
        raise ValueError(
                'Does not recognize the item "{}". Valid sections are '.\
                        format(line) + ', '.join(cls.__valid_sections) + '.')
        return None, None

    @staticmethod
    def _parse_natoms(args, timestep, data, start):
        """ Must not use timestep in this function """
        if args:
            raise RuntimeError("Non empty args")
        return int(data[start]), 1

    _parse_timestep = _parse_natoms
    #@classmethod
    #def _parse_timestep(cls, *args, **kwargs):
    #    return cls._parse_natoms(*args, **kwargs)

    @classmethod
    def _parse_bounds(cls, args, timestep, data, start):
        for arg in args:
            assert(arg == "pp") 
        ndim = len(args)
        assert(len(data) >= start+ndim)
        return cls._parse_formatted_section(["lbound", "hbound"], 
                [float, float],
                data[start:start+ndim]), ndim

    def _parse_atoms(self, args, timestep, data, start):
        natoms = self.trj[timestep]["natoms"]
        assert(len(data) >= start+natoms)
        return self._parse_formatted_section(args, 
                [self.data_types[arg] for arg in args],
                data[start:start+natoms]), natoms

    @staticmethod
    def _parse_formatted_section(args, types, data):
        rtn = []
        for line in data:
            elem = dict()
            for attr, type_, val in zip(args, types, line.split()):
                elem[attr] = type_(val)
            rtn.append(elem)
        return rtn

    # could be overridden by derived class
    __valid_keys = {
            "TIMESTEP" : "timestep", 
            "NUMBER OF ATOMS": "natoms", 
            "BOX BOUNDS" : "bounds", 
            "ATOMS" : "atoms",
            }
    __valid_sections = __valid_keys.keys()
    __data_types = {
            "id" : int,
            "mol" : int,
            "type" : int,
            "q" : float,
            "x" : float,
            "y" : float,
            "z" : float,
            "vx" : float,
            "vy" : float,
            "vz" : float,
            "fx" : float,
            "fy" : float,
            "fz" : float,
            }

if __name__ == "__main__":
    foo = LMPtrj()
    bar = LMPtrj("test/test.lammpstrj")
    foo.parse("test/test.lammpstrj")
    print(foo.trj)
    assert(foo.trj == bar.trj)
