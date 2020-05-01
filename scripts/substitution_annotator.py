from optparse import OptionParser
import re
import numpy as np


def getOptions():
    parser = OptionParser()
    parser.add_option("--input","-i", dest = "infile", help = "Input filtered sam file",
                      metavar = "FILE", type = "string", default = "")

    parser.add_option("--output_pre","-o", dest = "outfile", help = "Output prefix for sub samfiles",
                      metavar = "FILE", type = "string", default = "")

    (options, args) = parser.parse_args()
    return options

def main():
    options = getOptions()
    sam_file = options.infile
    out_prefix = options.outfile


    header = []

    reads = {}

   f_log = open(str(out_prefix + "_SNP.txt"),"w+")
   f_log.write("\t".join(("read_name","TT","TC","TG","TA","CT","CC","CG","CA","GT","GC","GG","GA","AT","AC","AG","AA"))+"\n")

    with open(sam_file, 'r') as f:
       for l,line in enumerate(f):

          line = line.strip()

          if line.startswith('@'):
             header.append(line)

          else:




             fields = line.split("\t")
             cigar = fields[5]
             seq = fields[9]
             qual = fields[10]
             name = fields[0]
             flag = int(fields[1])
             chr = fields[2]
             pos = fields[3]


             md_field = [m for m in fields[11:len(fields)] if m.startswith("MD")]
             if len(md_field) == 0:

                continue
             md = md_field[0].split(":")[2]
             md_ops = [x for x in re.split("[^a-zA-Z]*",md) if x]
             md_counts = [int(x) for x in re.split("[^0-9]*",md) if x]
             cigar_ops = [x for x in re.split("[^a-zA-Z]*",cigar) if x]
             cigar_counts = [int(x) for x in re.split("[^0-9]*",cigar) if x]
             q_ind = 0
             r_ind = 0
             md_temp = md
             tT=TC=TG=TA=0
             AT=AC=AG=tA=0
             CT=tC=CG=CA=0
             GT=GC=tG=GA=0



             try:
                for i in range(len(cigar_ops)):
                   m = cigar_ops[i]
                   n = cigar_counts[i]

                   if m == "M":
                      md_ind = 0
                      while md_ind < n:
                         token = re.search("^[0-9]*",md_temp).group(0)
                         if token != '':

                            md_temp = re.sub("^[0-9]*",'',md_temp)
                            if (md_ind+int(token)) > n:
                               md_temp = str(md_ind + int(token) - n) + md_temp
                               q_ind+=(n - md_ind)
                               r_ind+=(n - md_ind)
                               md_ind = n
                            else:
                               q_ind+=int(token)
                               r_ind+=int(token)
                               md_ind+=int(token)

                         else:
                            token = re.search("^[ATCGNnatcg]*",md_temp).group(0)
                            md_temp = re.sub("^[ATCGNnatcg]*",'',md_temp)



                           if len(token) > 0:
                               for j in range(len(token)):

                                  q_ind+=1
                                  md_ind+=1
                                  r_ind+=1
                                  s = seq[q_ind-1]
                                  r = token[j]
                                  q = ord(qual[q_ind-1])-33

                                  chr_pos = chr + "_" + str(int(pos)+r_ind-1)

                                  if flag < 256:

                                     if r == "T":
                                        if s == "A":
                                           TA+=1
                                        elif s == "C":
                                           TC+=1
                                        elif s == "G":
                                           TG+=1
                                     elif r == "A":
                                        if s == "T":
                                           AT+=1
                                        elif s == "C":
                                           AC+=1
                                        elif s == "G":
                                           AG+=1
                                     elif r == "C":
                                        if s == "T":
                                           CT+=1
                                        elif s == "A":
                                           CA+=1
                                        elif s == "G":
                                           CG+=1
                                     elif r == "G":
                                        if s == "T":
                                           GT+=1
                                        elif s == "A":
                                           GA+=1
                                        elif s == "C":
                                           GC+=1





                   elif m == "D":
                      r_ind+=n
                      token = re.search("^\^[ATCGNnatcg]*",md_temp).group(0)
                      md_temp = md_temp[n+1:len(md_temp)]
                   elif m == "S":
                      q_ind+=n
                   elif m == "I":
                      q_ind+=n
                   elif m == "N":
                      r_ind+=n

                if flag < 256:
                   tT=seq.count("T")-AT-CT-GT+TA+TC+TG
                   tA=seq.count("A")-TA-CA-GA+AT+AC+AG
                   tC=seq.count("C")-TC-AC-GC+CA+CG+CT
                   tG=seq.count("G")-AG-CG-TG+GA+GT+GC


                SNP = [tT,TC,TG,TA,CT,tC,CG,CA,GT,GC,tG,GA,AT,AC,AG,tA]
                if name in reads and flag < 256:
                   reads[name][0]= np.ndarray.tolist(np.add(reads[name][0],SNP))
                   reads[name].append(line)
                elif flag < 256:
                   reads[name] = [SNP,line]

             except:
                print("an error occured at line # :",l)

    for key,value in reads.items():
        f_log.write("\t".join([key]+[str(m) for m in value[0]])+"\n")

    f.close()
    f_log.close()
    f0 = open(str(out_prefix + "_TC_0.sam"),"w")
    f0.write("\n".join(header)+"\n")
    f0.write("\n".join([m[1] for m in reads.values()])+"\n")
    f0.write("\n".join([m[2] for m in reads.values() if len(m) == 3]))
    f0.close()
    f1 = open(str(out_prefix + "_TC_1.sam"),"w")
    f1.write("\n".join(header)+"\n")
    f1.write("\n".join([m[1] for m in reads.values() if m[0][1] >= 2])+"\n")
    f1.write("\n".join([m[2] for m in reads.values() if len(m) == 3 and m[0][1] >= 2]))
    f1.close()
    f2 = open(str(out_prefix + "_TC_2.sam"),"w")
    f2.write("\n".join(header)+"\n")
    f2.write("\n".join([m[1] for m in reads.values() if m[0][1] >= 4])+"\n")
    f2.write("\n".join([m[2] for m in reads.values() if len(m) == 3 and m[0][1] >= 4]))
    f2.close()
    f3 = open(str(out_prefix + "_TC_3.sam"),"w")
    f3.write("\n".join(header)+"\n")
    f3.write("\n".join([m[1] for m in reads.values() if m[0][1] >= 6])+"\n")
    f3.write("\n".join([m[2] for m in reads.values() if len(m) == 3 and m[0][1] >= 6]))
    f3.close()



if __name__ == '__main__':

   main()














