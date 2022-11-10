t = [1.0]
for i in range(1000000):
    a = t[-1] - t[-1]**2 /3
    t.append(a)
    print(t[-1])