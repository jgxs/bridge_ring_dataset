from requests_html import HTMLSession
import psycopg2
import psycopg2.extras
from bs4 import BeautifulSoup
import time
import re

tdata = [[]]
session = HTMLSession()
pages = ["","?page=1","?page=2","?page=3","?page=4","?page=5"]
for page in pages:
    r = session.get(f"https://www.kinase-screen.mrc.ac.uk/kinase-inhibitors{page}")  # 打开url网页 比如 driver.get("http://www.baidu.com")
    r.html.render(timeout=120)
    html = r.text
    tbody = BeautifulSoup(html, "html.parser").find("tbody")
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
    with open("test.html","a") as out:
        out.write(html)

print(len(tdata))
#for item in tdata:
#    print(str(item)[1:-1])

