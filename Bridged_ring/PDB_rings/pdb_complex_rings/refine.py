with open("bridge_complex.csv") as infor:
    infos = infor.read()

for item in infos.split("]"):
    print(item)