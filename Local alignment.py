f=open('blosum.txt','r')
x=f.readlines()
f.close()
l=[]
h=[""]

for i in x:
   singleRow = i.split(" ") 
   for j in singleRow:
       if j !="" and j !="\n":
           h.append(j)
   l.append(h)
   h=[]
def Local_Aligment_Protein(first_seq_pro, second_seq_pro):
    len_first_seq = len(first_seq_pro)
    len_second_seq = len(second_seq_pro)
    max_value = 0
    maxcell = (0, 0)
    Aligment = [[0 for first_seq_pro in range(len_first_seq + 1)] for second_seq_pro in range(len_second_seq + 1)]
    for i in range(len_second_seq + 1):
        Aligment[i][0] = 0
    for j in range(len_first_seq + 1):
        Aligment[0][j] = 0
    for i in range(1, len_second_seq + 1, 1):
        for j in range(1, len_first_seq + 1, 1):
           
                First_aligment_index = l[0].index(first_seq_pro[j-1])
                Second_aligment_index=l[0].index(second_seq_pro[i-1])
                score =l[First_aligment_index][Second_aligment_index] 
                Aligment[i][j] = max(Aligment[i - 1][j] - 10, Aligment[i][j - 1] - 10, Aligment[i - 1][j - 1] + int(score), 0)
                if max_value <= Aligment[i][j]:
                    max_value = Aligment[i][j]
                    maxcell = (i, j)

    print("This is Maximum value is", str(max_value))

    print("This is the Poistion of Maximum value in", str(maxcell))

    Gap_First_seq = ""
    Gap_Second_seq = ""
    Match = ""
    i, j = maxcell
    while (Aligment[i][j] > 0):
        up = Aligment[i - 1][j] - 10
        left = Aligment[i][j - 1] - 10
        
        First_aligment_index = l[0].index(first_seq_pro[j-1])
        Second_aligment_index=l[0].index(second_seq_pro[i-1])
            
        value=l[First_aligment_index][Second_aligment_index] 
        print("(",str(first_seq_pro[j-1]),",",str(second_seq_pro[i-1]),")")
            
           
        print("the value ",str(value))
        score = l[First_aligment_index][Second_aligment_index] 
        
        diagonal = Aligment[i - 1][j - 1] +int(score)
        if Aligment[i][j] == diagonal:
            Gap_First_seq += first_seq_pro[j - 1]
            Gap_Second_seq += second_seq_pro[i - 1]
           
            if first_seq_pro[j-1]==second_seq_pro[i-1]:
                Match += "|"
            else:
                Match += " "
            i -= 1
            j -= 1
        elif Aligment[i][j] == up:
            Gap_First_seq += "-"
            Match += " "
            Gap_Second_seq += second_seq_pro[i - 1]
            i -= 1
        else:
            Gap_First_seq += first_seq_pro[j - 1]
            Gap_Second_seq += "-"
            Match += " "
            j -= 1

    

    Gap_First_seq = Gap_First_seq[::-1]
    Match = Match[::-1]
    Gap_Second_seq = Gap_Second_seq[::-1]

    print (Gap_First_seq)
    print(Match)
    print(Gap_Second_seq)
   #----------------------this DNA Alignment function-------------
def Local_Aligment(first_seq, second_seq):
    len_first_seq=len(first_seq)
    len_second_seq=len(second_seq)
    max_value=0
    maxcell=(0,0)
    Aligment= [[0 for first_seq in range(len_first_seq+1)] for second_seq in range(len_second_seq+1)] 
    for i in range (len_second_seq+1):
        Aligment[i][0]=0
    for j in range (len_first_seq+1):
        Aligment[0][j] =0
    for i in range(1,len_second_seq+1,1):
        for j in range(1,len_first_seq+1,1):
            if first_seq[j-1] == second_seq[i-1]:
                score = 1
            else :
                score =-2
            Aligment[i][j] = max(Aligment[i-1][j]-1, Aligment[i][j-1]-1,Aligment[i-1][j-1]+score,0)
            if max_value<=Aligment[i][j]:
                 max_value=Aligment[i][j]
                 maxcell=(i,j)
                 
    print("This is Maximum value in Matrix is :",str(max_value))    
  
    print("This is the Poistion of Maximum value in :",str(maxcell))
    

    Gap_First_seq =""
    Gap_Second_seq=""
    Match=""
    i , j= maxcell
    while (Aligment[i][j] > 0):
        up = Aligment[i-1][j]-1
        left=Aligment[i][j-1]-1
        if  first_seq[j-1]==second_seq[i-1]:
            score=1
        else :
            score=-2
        diagonal = Aligment[i-1][j-1]+score
        if Aligment[i][j]==diagonal:
            Gap_First_seq +=first_seq[j-1]
            Gap_Second_seq+=second_seq[i-1]
            if score==1:
                Match+="|"
            else:
                Match+=" "
            i-=1
            j-=1
        elif Aligment[i][j] == up:
            Gap_First_seq +="-"
            Match+=" "
            Gap_Second_seq+=second_seq[i-1]
            i-=1
        else:
            Gap_First_seq +=first_seq[j-1]
            Gap_Second_seq+="-"
            Match+=" "
            j-=1

   
    Gap_First_seq = Gap_First_seq[::-1]
    Match=Match[::-1]
    Gap_Second_seq=Gap_Second_seq[::-1]
    
    print (Gap_First_seq )
    print(Match)
    print(Gap_Second_seq)
    
 
# calling the function
choice=""
first_seq=""
second_seq =""
first_seq_pro = ""
second_seq_pro = ""
print("Your choice may be DNA ot Protein")
choice=input("Enter choice")
if choice=="DNA" or choice=="dna":
   print("Enter the First Sequence of DNA") 
   first_seq=input()
   print("Enter the Second Sequence of DNA") 
   second_seq=input()
   Local_Aligment(first_seq, second_seq)
elif choice=="Protein" or choice=="protein": 
    print("Enter the First Sequence of Protein") 
    first_seq_pro=input()
    print("Enter the  second Sequence of Protein") 
    second_seq_pro=input()
    Local_Aligment_Protein(first_seq_pro,second_seq_pro)
    
    
else:
 print("Invalid choice")

