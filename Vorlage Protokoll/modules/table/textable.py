import uncertainties
def latex_tab(data=[[1,2,3],[42,42,42]],names=["col1","col2"],filename="test.tex",caption="Caption",label="test", dec_points=[0,2]):
    try:
        for i in range(len(data)-1):
            if not(len(data[i])==len(data[i+1])):
                raise TypeError("Arrays must have same dimension!"+" Error occured in comparison of data["+str(i)+"]"+"and data["+str(i+1)+"]")



        if not(len(names)==len(data) and len(dec_points)==len(data)):
            raise TypeError("data and names and dec_points must have same dimension! "+"len(data)= "+str(len(data))+"; len(names)= "+str(len(names)) +"; len(dec_points)= "+str(len(dec_points)))
    except TypeError: raise
    else:
	texfile = open(filename,"w")
        texfile.write("\\begin{table}\n");
        texfile.write(" \\caption{"+caption+"}\n");
        texfile.write(" \\label{tab:"+label+"}\n")
        texfile.write(" \\centering\n")
        texfile.write("\\sisetup{table-format=1.1}")
        texfile.write(" \\begin{tabular}{")
        for col in data:
            texfile.write("S")
        texfile.write("}\n");
        texfile.write(" \\toprule \n    ");
        for i in range(len(names)-1):
            texfile.write("{"+names[i]+"}& ")
        texfile.write("{"+names[len(names)-1]+"}")
        texfile.write(" \\\\\n");
        texfile.write("     \\midrule\n");
        for i in range(len(data[0])):
            texfile.write("     ")
            for j in range(len(data)-1):

                if isinstance(data[j][i],uncertainties.core.Variable):
                    data[j][i]=(str(data[j][i])).replace('+/-',' \pm ')
                    texfile.write('$\\num{'+ data[j][i]+'}$')
                else:
                    texfile.write(("{:10.%df}"%dec_points[j]).format(data[j][i])+" & ")
            if isinstance(data[len(data)-1][i],uncertainties.core.Variable):
                data[len(data)-1][i]=(str(data[len(data)-1][i])).replace('+/-','\pm')
                texfile.write('$\\num{'+ data[len(data)-1][i]+'}$')
            else:
                texfile.write(("{:10.%df}"%dec_points[len(dec_points)-1]).format(data[len(data)-1][i]))
            texfile.write(" \\\\\n")
        texfile.write(" \\bottomrule\n");
        texfile.write(" \\end{tabular}\n");
        texfile.write("\\end{table}");
        texfile.close()




if __name__ == "__main__":
    arr1=[1,2,3,4,5]
    arr2=[0.8,12,234,42,1.28]
    import textable
    print("42 is the answer to life the universe and everything!")
    print("Running this module as script generates a sample table.tex")
    latex_tab([arr1,arr2],["col1","col2"],"table.tex",'caption','sample')

