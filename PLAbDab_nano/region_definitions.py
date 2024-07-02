import numpy as np

# All are north region definitions in imgt numbering
reg_def = dict()
reg_def["fwH"] = [6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                50, 51, 52, 53, 54, 55, 66, 67, 68, 69, 70, 71, 72, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87,
                88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 118, 119, 120, 121, 122, 123]

reg_def["CDRH1"] = [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38]
reg_def["CDRH2"] = [56, 57, 58, 59, 60, 61, 62, 63, 64, 65]
reg_def["CDRH3"] = [105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117]

reg_def["CDR_H_all"] = reg_def["CDRH1"] + reg_def["CDRH2"] + reg_def["CDRH3"]
reg_def["CDR_all"] = reg_def["CDR_H_all"] 
reg_def["fw_all"] = reg_def["fwH"] 

# numba does not like lists very much
reg_def = {x: np.array(reg_def[x]) for x in reg_def}

# Store these to use in jit
reg_def_CDR_all = reg_def["CDR_all"]
reg_def_fw_all = reg_def["fw_all"]