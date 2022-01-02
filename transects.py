file = open("transects.txt", "w")
file.write("# line X Y\n")
file.write("1 -100 -100\n")
X = -100

for i in range(0, 11):
    Y = -100 + 20 * i
    X = -X
    
    line = "1 {0} {1}\n".format(X, Y)
    file.write(line)
    
    if i != 10:
        line = "1 {0} {1}\n".format(X, Y + 20)
        file.write(line)

file.close()
