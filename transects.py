file = open("transects.txt", "w")
file.write("# line X Y\n")

for i in range(0, 11):
    Y = -100 + 20 * i
    
    if i % 2 == 0:
        X = -100    
    elif i % 2 != 0:
        X = 100
        
    line = "1 {0} {1}\n".format(X, Y)
    file.write(line)
    line = "1 {0} {1}\n".format(-X, Y)
    file.write(line)

    if i != 10:
        line = "1 {0} {1}\n".format(-X, Y + 20)
        file.write(line)

file.close()
