from configparser import ConfigParser
import os

cnf = ConfigParser()
dirname = os.path.dirname(__file__)
cnf.read(os.path.join(dirname, "./config.ini"))

ENTRAZ_EMAIL = cnf.get("BIOPYTHON", "EMAIL", fallback="me21m014@technikum-wien.at")
GENE_NAME = cnf.get("BIOPYTHON", "GENE")
GENE_NAME_ALTERNATIVE_1 = cnf.get("BIOPYTHON", "GENE_ALTERNATIVE_1")
GENE_NAME_ALTERNATIVE_2 = cnf.get("BIOPYTHON", "GENE_ALTERNATIVE_2")
GENE_NAME_ALTERNATIVE_3 = cnf.get("BIOPYTHON", "GENE_ALTERNATIVE_3")
GENE_NAME_ALTERNATIVE_4 = cnf.get("BIOPYTHON", "GENE_ALTERNATIVE_4")