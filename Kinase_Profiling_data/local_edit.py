from requests_html import HTMLSession
from bs4 import BeautifulSoup
import psycopg2
import psycopg2.extras

def getsmiles(compound_name):
    with psycopg2.connect(dbname="chembl_29", user="user", password="user", host="192.168.54.19") as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cursor:
            cursor.execute("select molregno from compound_records where compound_name=%s;",(compound_name,))
            try:
                Molrengno = cursor.fetchall()[0][0]
            except:
                return False
            cursor.execute("select chembl_id from MOLECULE_DICTIONARY where molregno=%s;",(Molrengno,))
            Chembl_id = cursor.fetchall()[0][0]
            cursor.execute("select canonical_smiles  from COMPOUND_STRUCTURES where molregno=%s;",(Molrengno,))
            try:
                SMILES = cursor.fetchall()[0][0]
            except:
                return False
        return [Chembl_id,SMILES]

tdata = [[]]
with open("/home/chengyj/kinase_work/data/test.html") as r:
 # 打开url网页 比如 driver.get("http://www.baidu.com")
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


for item in tdata[0:-1]:
    id_smi = getsmiles(item[0])
    if id_smi: 
        print(f"{item[0]:<60} {item[4]:<6} {item[3]:<20} {id_smi[0]:<15} {id_smi[1]}")
    else:
        try:
            print(f"{item[0]:<60} {item[4]:<6} {item[3]:<20} https://www.kinase-screen.mrc.ac.uk/kinase-inhibitors/{item[2]}")
        except:
            print(item)

