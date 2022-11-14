bch = ""
with open("menu2.csv") as aaa:
    for line in aaa:
        if "bch" in line:
            bch = line
        else:
            if line[0:8] == bch[0:8]:
                pass
            else:
                print(line[0:-1])