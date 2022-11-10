from urllib.parse import quote
import requests

smi = {}
with open("/home/chengyj/kinase_work/Bridged_ring/chembl_rings/KI_chembl/spiro_bridged_ki.csv") as smi_file:
    for line in smi_file:
        smi[line.split()[0]] = line.split()[1]
with open("/home/chengyj/kinase_work/Bridged_ring/chembl_rings/KI_chembl/data_refined/spiro_ki.csv") as smi_file:
    for line in smi_file:
        smi[line.split()[0]] = line.split()[1]
with open("/home/chengyj/kinase_work/Bridged_ring/chembl_rings/KI_chembl/data_refined/bridged_ki.csv") as smi_file:
    for line in smi_file:
        smi[line.split()[0]] = line.split()[1]

url = """{"query":{"type":"terminal","service":"chemical","parameters":{"value":"SMILES_INPUT","type":"descriptor","descriptor_type":"SMILES","match_type":"graph-strict"}},"return_type":"entry"}"""
for key in smi:
    print(key)
    url_re = "https://search.rcsb.org/rcsbsearch/v2/query?json="+url.replace("\"",quote("\"")).replace("{",quote("{")).replace("}",quote("}")).replace(":",quote(":")).replace(",",quote(",")).replace("SMILES_INPUT",smi[key])
    r = requests.get(url_re) 
    if r.status_code == 200:
        result = [item.replace("\"","").replace("\n","").replace(" ","").replace("{","").replace("}","").replace("]","").replace(",score","").replace(":","",1) for item in r.text.split("identifier")[1:]]
        with open("/home/chengyj/kinase_work/Bridged_ring/pdb_rings/bridge_complex.csv","a") as output:
            output.write(f"{key:>8} {result}")