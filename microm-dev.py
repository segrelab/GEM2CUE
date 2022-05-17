import micom
from micom.data import test_taxonomy

taxonomies = [test_taxonomy(n=n) for n in range(2, 12)]

print(taxonomies[2])