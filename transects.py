file = open("transects.txt", "w")
file.write("# line X Y\n")
X = 100

for i in range(1, 12):
    Y = -120 + 20 * i
    X = -X
    
    line = "{0} {1} {2}\n".format(i, X, Y)
    file.write(line)
    line = "{0} {1} {2}\n".format(i, -X, Y)
    file.write(line)
    
file.close()
