def latex_tab(data=[[1,2,3],[42,42,42]],names=["col1","col2"],filename="test.tex",caption="Caption",label="test"):
    try:
        for i in range(len(data)-1):
            if not(len(data[i])==len(data[i+1])):
                raise TypeError("Arrays must have same dimension!"+" Error occured in comparison of data["+str(i)+"]"+"and data["+str(i+1)+"]")



        if not(len(names)==len(data)):
            raise TypeError("data and names must have same dimension! "+"len(data)= "+str(len(data))+"; len(names)= "+str(len(names)))
    except TypeError: raise
    else:
        texfile = open(filename,"w")
        texfile.write("\\begin{table}\n");
        texfile.write(" \\caption{"+caption+"}\n");
        texfile.write(" \\label{tab:"+label+"}\n")
        texfile.write(" \\centering\n")
        texfile.write(" \\begin{tabular}{")
        for col in data:
            texfile.write("c")
        texfile.write("}\n");
        texfile.write(" \\toprule \n    ");
        for i in range(len(names)-1):
            texfile.write(names[i]+" & ")
        texfile.write(names[len(names)-1])
        texfile.write(" \\\\\n");
        texfile.write("     \\midrule\n");
        for i in range(len(data[0])):
            texfile.write("     ")
            for j in range(len(data)-1):
                texfile.write(str(data[j][i])+" & ")
            texfile.write(str(data[len(data)-1][i]))
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
