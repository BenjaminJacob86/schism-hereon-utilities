

def techit(fname='test.tex',title='manuscript',text=''):

	from os import listdir
	from subprocess import call


	file=open(fname,'w')

	file.write("\\documentclass[a4paper]{article}\n")
	file.write("\\usepackage[a4paper]{geometry}\n")
	file.write("\\parindent0pt \n \\textwidth16cm \n \\textheight22cm \n")

	#file.write("\\usepackage{graphicx}\n")
	file.write("\\usepackage[pdftex]{graphicx}\n")
	file.write("\\usepackage{epsfig}\n")
	file.write("\\usepackage{epstopdf}\n")


	file.write("\\begin{document}\n")

	file.write("\\title{"+title+"}\n")
	file.write("\\maketitle\n")
	
	file.write("\n"+text+"\n")

	count=1
	files=sorted(listdir('./'))
	for picture in files:
		if ".png" in picture or ".eps" in picture:
			count=count+1
			print(picture)
			file.write("\n")
			file.write("\\begin{figure} \n \\centering \n")
			file.write(" \\includegraphics[width=1.00\\textwidth]{%s} \n" %(picture))
			file.write("\\end{figure}\n")
			if count%2==0:
				file.write("\\clearpage\n")

	file.write("\\end{document}")

if __name__ == '__main__':
	techit(fname='manuscript.tex',title='manuscript',text='')