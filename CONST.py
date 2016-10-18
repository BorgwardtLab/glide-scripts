# CONST.py -- Provide the dictionary necessary to the conversion between
# tile IDs and indices

# blockToIndex is computed with compute_block_to_index.py
# using the appropriate TILESIZE.
# Its purpose is to map the file tiles (created by splitting the input files)
# back to the corresponding SNP indices

# This dictionary is used by compute_all_pvalues.py

# contact: chloe-agathe.azencott@mines-paristech.fr



# TILESIZE = 24576
blockToIndex = {'aa': 0, 'ac': 49152, 'ab': 24576, 'ae': 98304, 'ad': 73728, 'ag': 147456, 'af': 122880, 'ai': 196608, 'ah': 172032, 'ak': 245760, 'aj': 221184, 'am': 294912, 'al': 270336, 'ao': 344064, 'an': 319488, 'aq': 393216, 'ap': 368640, 'as': 442368, 'ar': 417792, 'au': 491520, 'at': 466944, 'aw': 540672, 'av': 516096, 'ay': 589824, 'ax': 565248, 'az': 614400, 'bd': 712704, 'be': 737280, 'bf': 761856, 'bg': 786432, 'ba': 638976, 'bb': 663552, 'bc': 688128, 'bl': 909312, 'bm': 933888, 'bn': 958464, 'bo': 983040, 'bh': 811008, 'bi': 835584, 'bj': 860160, 'bk': 884736, 'bt': 1105920, 'bu': 1130496, 'bv': 1155072, 'bw': 1179648, 'bp': 1007616, 'bq': 1032192, 'br': 1056768, 'bs': 1081344, 'bx': 1204224, 'by': 1228800, 'bz': 1253376}


# # TILESIZE = 10000
# blockToIndex = {'aa': 0, 'ac': 20000, 'ab': 10000, 'ae': 40000, 'ad': 30000, 'ag': 60000, 'af': 50000, 'ai': 80000, 'ah': 70000, 'ak': 100000, 'aj': 90000, 'am': 120000, 'al': 110000, 'ao': 140000, 'an': 130000, 'aq': 160000, 'ap': 150000, 'as': 180000, 'ar': 170000, 'au': 200000, 'at': 190000, 'aw': 220000, 'av': 210000, 'ay': 240000, 'ax': 230000, 'az': 250000, 'bd': 290000, 'be': 300000, 'bf': 310000, 'bg': 320000, 'ba': 260000, 'bb': 270000, 'bc': 280000, 'bl': 370000, 'bm': 380000, 'bn': 390000, 'bo': 400000, 'bh': 330000, 'bi': 340000, 'bj': 350000, 'bk': 360000, 'bt': 450000, 'bu': 460000, 'bv': 470000, 'bw': 480000, 'bp': 410000, 'bq': 420000, 'br': 430000, 'bs': 440000, 'bx': 490000, 'by': 500000, 'bz': 510000}
