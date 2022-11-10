import psycopg2
import psycopg2.extras
from bs4 import BeautifulSoup

class ligand_info:
    def __init__(self,line):
        self.name = line[0:60].rstrip().lstrip()
        self.info = line[60:].split()
        self.weight = float(self.info[0])
        self.formula = self.info[1]
        if "CHEMBL" in line:
            self.chemblid = line[line.find("CHEMBL"):].split()[0]
            self.smiles = "Not found" 
            self.getsmile()
        else:
            self.smiles = self.info[-2]
    def getsmile(self):
        with psycopg2.connect(dbname="chembl_29", user="user", password="user", host="192.168.54.19") as conn:
            with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
                cursor.execute("select molregno from MOLECULE_DICTIONARY where molecule_dictionary.chembl_id=%s;",(self.chemblid,))
                self.Molrengno = cursor.fetchall()[0][0]
                cursor.execute("select canonical_smiles  from COMPOUND_STRUCTURES where molregno=%s;",(self.Molrengno,))
                self.smiles = cursor.fetchall()[0][0]

datas = {}
with open("with") as all_data:
    for line in all_data:
        ligand = ligand_info(line)
        datas[ligand.name] = ligand

tdata = [[]]
with open("/home/chengyj/kinase_work/data/test.html") as r:
    html = r
    tbodys = BeautifulSoup(html, "html.parser").find_all("tbody")
    for tbody in tbodys:
        trows = tbody.find_all("tr")
        for trow in trows:
            for item in trow.find_all("td"):
                text = item.text.rstrip().lstrip()
                if text != "":
                    tdata[-1].append(text)
                if item.a:
                    if "screening-compounds" in item.a.get("href") or "jpg" in item.a.get("href"): 
                        href = item.a.get("href")
                        tdata[-1].append(href)
            tdata.append([])

for item in tdata:
    try:
        print(f"{datas[item[0]].smiles:<180} {item[0]:<60} {item[2].split('/')[-1]}")
    except:
        pass

