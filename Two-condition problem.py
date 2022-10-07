#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('display', 'latex')
import sympy as sp
from itertools import combinations
import re


# In[2]:


import string
def arb_vector_field( Manifold, chart, frame, Name=None):
    if dim(Manifold)>26:
        return ERROR
    
    alphabet_string = string.ascii_lowercase
    alphabet_list = list(alphabet_string)
    
    vectorlist = [None]*dim(Manifold)
    
    for i in range(dim(Manifold)):
        vectorlist[i] = alphabet_list[i]+"_"+Name
        
    #func = [None]*dim(M)
    for i in range (len(vectorlist)):
        locals()[vectorlist[i]] = var(vectorlist[i])
    for i in range (len(vectorlist)):
        vectorlist[i]=locals()[vectorlist[i]]
    
    if Name==None:
        return (Manifold.vector_field({frame:vectorlist}))
    else:
        return (Manifold.vector_field({frame:vectorlist},name = Name),vectorlist)


# In[3]:


def Torsion_checker(Ideal, Independence_Conditions, vector_field_1, vector_field_2, Manifold, variables, cond):
    IC=[]
    dA1 = Ideal[3]
    dA2 = Ideal[4]
    dA3 = Ideal[5]
    constraints = [None]*8
    for k in range(len(variables[2:6])):
        constraints[k] = (variables[k+2]!=0)
    for l in range(len(variables[8:12])):
        constraints[l+len(variables[2:6])] = (variables[l+8]!=0) 
    
    for j in range(len(Independence_Conditions)):
        IC.append(vector_field_1(Independence_Conditions[0]).expr()-1)
        IC.append(vector_field_2(Independence_Conditions[0]).expr())
    
    IC.append(vector_field_1(Independence_Conditions[1]).expr())
    IC.append(vector_field_2(Independence_Conditions[1]).expr()-1)
    
    for k in range(len(Ideal)/2):
        IC.append(vector_field_1(Ideal[k]).expr())
        IC.append(vector_field_2(Ideal[k]).expr())
    
    #soln1 = solve(IC,variables)[0]
    
    for l in range(len(Ideal)/2,len(Ideal)):
        IC.append((Ideal[l])(vector_field_1,vector_field_2).expr())
    
    soln2 = solve(IC+cond,tuple(i for i in variables))[0]
    if soln2 == []:
        print("There exists torsion in the system")
    else:
        print("The torsion in the system is absorbable")
    con = 1
    lst_of_cond = []
    for i in range(len(soln2)):
        if bool(soln2[i].right_hand_side().denominator()==1)==False : 
            if soln2[i].right_hand_side().denominator() not in lst_of_cond:
                lst_of_cond.append(soln2[i].right_hand_side().denominator())
                lst_of_cond.append(soln2[i].right_hand_side().numerator())
    
                print("{0} Dimensional Integral element exists iff {1} = 0, {2} = 0".format(con,(soln2[i].right_hand_side().denominator()),
                                                                                                     soln2[i].right_hand_side().numerator()))
                con+=1
    for i in range(len(soln2)):
        soln2[i] = soln2[i].right_hand_side()                     
    return soln2[0:6],soln2[6:12], lst_of_cond
    #return(soln2)


# In[4]:


def Cartan_character_Calculator( M, Ideal, X, Y):
    weg = wedge_max(M,Ideal)
    
    TableuaxX = []
    TableuaxY = []
    tallyX = 0
    for diff in Ideal:
        tally = 0
        if diff.degree()==2:
            diff_form = X.contract(diff)
            for i in range(len(weg)):
                tmp = diff_form.wedge(weg[i])
                if tmp!=0:
                    tally += 1
            if tally == len(weg):
                tallyX += 1
                TableuaxX.append(diff_form)
    tallyY = 0            
    for diff in Ideal:
        tally = 0
        if diff.degree()==2:
            diff_form = Y.contract(diff)
            for i in range(len(weg)):
                tmp = diff_form.wedge(weg[i])
                if tmp!=0:
                    tally += 1
            if tally == len(weg):
                tallyY +=1
                TableuaxY.append(diff_form)
    
    Tableuax = TableuaxX + TableuaxY
    if len(Tableuax)==0:
        print("s_1=0")
    elif len(Tableuax)==1:
        print("s_1=1")
    else:
        for i in range(1,len(TableuaxX)):           
            for j in range(0,i-1):
                if TableuaxX[j].wedge(TableuaxX[i])==0:
                    TableuaxX[i]=1
                    break
        for k in range(len(TableuaxX)):
            if TableuaxX[k]==1:
                TableuaxX.remove(TableuaxX[k])
                
        for i in range(1,len(TableuaxY)):           
            for j in range(0,i-1):
                if TableuaxY[j].wedge(TableuaxY[i])==0:
                    TableuaxY[i]=1
                    break
        for k in range(len(TableuaxY)):
            if TableuaxY[k]==1:
                TableuaxY.remove(TableuaxY[k])
                
        if len(TableuaxX)>=len(TableuaxY):
            print("s_1 = {0}".format(len(TableuaxX)))
            s_2 = len(TableuaxY)
            for i in range(len(TableuaxY)):
                for j in range(len(TableuaxX)):
                    tmp = TableuaxY[i].wedge(TableuaxX[j])
                    if tmp == 0:
                        s_2 = s_2-2
                        break
                        
        elif len(TableuaxY)<len(TableuaxY):
            print("s_1 = {0}".format(len(TableuaxY)))
            s_2 = len(TableuaxX)
            for i in range(len(TableuaxX)):
                for j in range(len(TableuaxY)):
                    tmp = TableuaxX[i].wedge(TableuaxY[j])
                    if tmp == 0:
                        s_2 = s_2-1
                        break
                        
        print("s_2= {0}".format(s_2))


# In[5]:


def wedge_max(M,Ideal):
    weg = []
    weg.append(Ideal[0])
    for arg in range(1,len(Ideal)):
        for k in range(0,len(weg)):
            wedg = weg[k].wedge(Ideal[arg])
            if wedg!=0 and (wedg).degree()!=dim(M):
                weg[k] = wedg
                break
            elif k==len(weg)-1:
                weg.append(Ideal[arg])
                break 
    return (weg)


# In[6]:


def Independence_Conditions_Checker(M,list_of_indep, diff_form):
    for i in range(len(list_of_indep)):
        if diff_form.wedge(list_of_indep[i])==0:
            return("Contracted Differential Form is not independent of the list of differential forms")
    list_of_cond = []
    res=[]
    lst=[]
    for i in range(0,dim(M)):
            lst.append(i)
    for i in range(0,len(list_of_indep)):
        
        wedg = list_of_indep[i].wedge(diff_form)    
        degg = list_of_indep[i].wedge(diff_form).degree()

        comb = list(combinations(lst, int(degg)))
        
        for j in comb:
            if list_of_ind[i].wedge(diff_form)[j]!=0:
                if list_of_ind[i].wedge(diff_form)[j].expr().variables() == ():
                    return(LatexExpr("Differential \: form \; {} \: is \: always \: independent \: of \: differential \: forms \: given".format(latex(diff_form.parent()))))
                if list_of_ind[i].wedge(diff_form)[j] not in res and -list_of_ind[i].wedge(diff_form)[j] not in res:
                    res.append(list_of_ind[i].wedge(diff_form)[j])
            
        
    return (res)


# In[7]:


def loop_rec(y, degg, n,wedg,res ):
    inp = [0]*n
    if n >= 1:
        for x in range(y, degg-1):
            inp[x] = y
            loop_rec(y+1, n - 1, degg,wedg,res)
    else:
        if wedg[inp]==1 or wedg[inp]==-1:
            return("always")

        if wedg[inp] not in res and -(wedg[inp]) not in res:
            return(res.append(wedg[inp]))
    return(res)


# In[8]:


R6= Manifold(6,"R^{6}")
chartR6.<x,y,u,u_x,u_y,u_yy> = R6.chart()

eA = chartR6.frame()
 
F = function('F',nargs=6)(*chartR6)
G = function('G',nargs=6)(*chartR6)

A4 = R6.one_form({eA:[1,0,0,0,0,0]}, name='\Omega_4')
A5 = R6.one_form({eA:[0,1,0,0,0,0]}, name='\Omega_4')

A1 = R6.one_form({eA:[-u_x,-u_y,1,0,0,0]}, name='\Omega_1')
A2 = R6.one_form({eA:[-F,-G,0,1,0,0]}, name='\Omega_2')
A3 = R6.one_form({eA:[-G,-u_yy,0,0,1,0]}, name='\Omega_3')

X, varx  = arb_vector_field( R6, chartR6, eA, "X")

Y,varY = arb_vector_field( R6, chartR6, eA, "Y")



dA2=A2.exterior_derivative()
dA1=A1.exterior_derivative()
dA3=A3.exterior_derivative()
Ideal = [A1,A2,A3,dA1,dA2,dA3]
Independence = [A4,A5]


# In[9]:


X, varx  = arb_vector_field( R6, chartR6, eA, "X")
Y,varY = arb_vector_field( R6, chartR6, eA, "Y")


# In[10]:


X[:], Y[:], cond = Torsion_checker(Ideal, Independence, X, Y, R6, varx+varY, [])


# In[11]:


dA2(X,Y).display()


# In[12]:


cond[0]


# In[13]:


dA2.display()


# In[14]:


dA2[0,5] = diff(G,u_yy)


# In[15]:


X, varx  = arb_vector_field( R6, chartR6, eA, "X")
Y,varY = arb_vector_field( R6, chartR6, eA, "Y")
X[:], Y[:], cond = Torsion_checker(Ideal, Independence, X, Y, R6, varx+varY, [])


# In[16]:


X[5], Y[5]


# In[17]:


# dG/du_yy = 1
dA2[0,5] = 1
dA2[1,5] = 1
dA3[0,5] = 1
dA3.display()


# In[18]:


X[5].display(), Y[5].display()


# In[19]:


locals()["a"] = var('a')
X[5] = a
Y[5] = a


# In[20]:


dA2(X,Y).display()


# In[21]:


#dG/dx
cond_dg_dx = -u_y*diff(F, u) - G*diff(F, u_x) - u_yy*diff(F, u_y) + u_x*diff(G, u) + F*diff(G, u_x) + G*diff(G, u_y)
dA2[0,1] = cond_dg_dx
dA2.display()


# In[22]:


dA2(X,Y).display()


# In[23]:


dA3(X,Y).display()


# In[24]:


cond_dg_dy = -u_y*diff(G,u)-G*diff(G,u_x)- u_yy*diff(G, u_y)
dA3[0,1] = cond_dg_dy


# In[25]:


dA3.display()


# In[26]:


dA3(X,Y).display()


# In[ ]:


#The torsion is now absorbable with one free variable on the vector fields


# In[28]:


X.display(), Y.display()


# In[ ]:


#The system is involution


# In[29]:


Ideal = [A1, A2, A3, dA1, dA2, dA3]
Cartan_character_Calculator( R6, Ideal, X, Y)

