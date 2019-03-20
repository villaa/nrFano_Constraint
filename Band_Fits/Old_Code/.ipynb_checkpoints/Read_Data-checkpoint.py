def Read_File(file_name):
    with open(file_name, 'r') as data:
        Enr = []
        Yield = []
        for line in data:
            p = line.split()
            Enr.append(float(p[0]))
            Yield.append(float(p[1]))

    return Enr, Yield



