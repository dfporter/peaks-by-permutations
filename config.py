import ConfigParser
#import sys
import os

def config():
    Config = ConfigParser.ConfigParser()
    src_dir = os.path.dirname(os.path.realpath(__file__))
    Config.read(src_dir + "/config.ini")
    lib = ConfigSectionMap('library', Config)
    return lib

 
def ConfigSectionMap(section, Config):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1
