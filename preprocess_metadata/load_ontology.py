from onto_lib_py3 import load_ontology as lo

ONT_NAME_TO_ONT_ID = {
    "EFO_CL_DOID_UBERON_CVCL": "17",
    "UO": "7"    
}
ONT_ID_TO_OG = {x: lo.load(x)[0] for x in ONT_NAME_TO_ONT_ID.values()}

def the_ontology():
    return ONT_ID_TO_OG['17']

def unit_ontology():
    return ONT_ID_TO_OG['7']
