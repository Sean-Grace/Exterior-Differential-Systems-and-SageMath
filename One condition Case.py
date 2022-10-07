#!/usr/bin/env python
# coding: utf-8

# In[5]:


get_ipython().run_line_magic('display', 'latex')
import sympy as sp
from itertools import combinations


# In[6]:


import string
def arb_vector_field( Manifold, chart, frame, Name=None):
    if dim(Manifold)>26:
        return ERROR
    
    alphabet_string = string.ascii_lowercase
    alphabet_list = list(alphabet_string)
    
    vectorlist = [None]*dim(Manifold)
    
    if Name==None:
        vectorlist[i] = alphabet_list[i]
    else:
        for i in range(dim(Manifold)):
            vectorlist[i] = alphabet_list[i]+"_"+Name
        
    #func = [None]*dim(M)
    for i in range (len(vectorlist)):
        locals()[vectorlist[i]] = var(vectorlist[i])
    for i in range (len(vectorlist)):
        vectorlist[i]=locals()[vectorlist[i]]
    
    if Name==None:
        return (Manifold.vector_field({frame:vectorlist}), vectorlist)
    else:
        return (Manifold.vector_field({frame:vectorlist},name = Name), vectorlist)


# In[7]:


def Condition_on_arb_vector_field(Manifold, frame, X, Omega, result, var_list):
    Eval = Omega(X).expr()
    Xl = X[:]
    rem_var = []
    for i in range(len(var_list)):
        soln = sp.solveset(Eval-result, var_list[i])
        for num in soln:
            Xl[i] = num
    return (Manifold.vector_field({frame:Xl}, name = my_name(X)))


# In[8]:


def Condition_on_2_arb_vectorfield(Manifold , frame, X, Y, two_form, result, var_list):
    Eval = two_form(X,Y).expr()
    Xl = X[:]
    Yl = Y[:]
    for i in range(len(var_list)):
        soln = sp.solveset(Eval-result,var_list[i])
        for num in soln:
            Xl[i]=num
    return( Manifold.vector_field({frame: Xl}, my_name(X)))


# In[9]:


def wedge_max(M, Ideal):
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


# In[10]:


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


# In[11]:


def Independence_Conditions_Checker( M, list_of_indep, diff_form, v_field):
    for i in range(len(list_of_indep)):
        if diff_form.wedge(list_of_indep[i])==0:
            return("Differential Form is not independent of the list of differential forms")
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

            if list_of_ind[i].wedge(diff_form)[j].expr().variables() == ():
                print(LatexExpr("Differential \: form \; {} \: is \: always \: independent \: of \: differential \: forms \: given".format(latex(diff_form.parent()))))
                return(list_of_indep.append(diff_form))
            if list_of_ind[i].wedge(diff_form)[j] not in res and -list_of_ind[i].wedge(diff_form)[j] not in res:
                res.append(list_of_ind[i].wedge(diff_form)[j])
            
        
    return (res)


# In[12]:


def my_name(p):
    for name, value in globals().items():
        if id(value) == id(p):
            return name


# In[13]:


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
        print(0)
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
                        s_2 = s_2-1
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


# In[14]:


R7= Manifold(7,"R^{7}")
chartR7.<x,y,u,u_x,u_y,u_xy,u_yy> = R7.chart()

eA = chartR7.frame()  
    
F = function('F',nargs=7)(*chartR7)

A4 = R7.one_form({eA:[1,0,0,0,0,0,0]}, name='\Omega_4')
A5 = R7.one_form({eA:[0,1,0,0,0,0,0]}, name='\Omega_5')

A1 = R7.one_form({eA:[-u_x,-u_y,1,0,0,0,0]}, name='\Omega_1')
A2 = R7.one_form({eA:[F,-u_xy,0,1,0,0,0]}, name='\Omega_2')
A3 = R7.one_form({eA:[-u_xy,-u_yy,0,0,1,0,0]}, name='\Omega_3')

X, varx = arb_vector_field( R7, chartR7, eA, "X")

Y, vary = arb_vector_field( R7, chartR7, eA, "Y")



dA2=A2.exterior_derivative()
dA1=A1.exterior_derivative()
dA3=A3.exterior_derivative()


# In[15]:


X.display(), Y.display()


# In[16]:


X =Condition_on_arb_vector_field(R7,eA,X,A4,1,varx)
X =Condition_on_arb_vector_field(R7,eA,X,A5,0, varx)
X =Condition_on_arb_vector_field(R7,eA,X,A1,0,varx)
X =Condition_on_arb_vector_field(R7,eA,X,A2,0,varx)
X =Condition_on_arb_vector_field(R7,eA,X,A3,0,varx)


# In[17]:


X.display()


# In[18]:


Y =Condition_on_arb_vector_field(R7,eA,Y,A4,0,vary)
Y =Condition_on_arb_vector_field(R7,eA,Y,A5,1,vary)
Y =Condition_on_arb_vector_field(R7,eA,Y,A1,0,vary)
Y =Condition_on_arb_vector_field(R7,eA,Y,A2,0,vary)
Y =Condition_on_arb_vector_field(R7,eA,Y,A3,0,vary)


# In[19]:


Y.display()


# In[20]:


X=Condition_on_2_arb_vectorfield(R7,eA,X,Y,dA3,0,varx)
X=Condition_on_2_arb_vectorfield(R7,eA,X,Y,dA2,0,varx)


# In[21]:


X.display()


# In[22]:


Y.display()


# In[23]:


Cartan_character_Calculator(R7, [A1,A2,A3,dA1,dA2,dA3] , X, Y)

