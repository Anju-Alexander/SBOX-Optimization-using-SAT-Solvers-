# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License version 2 as published
# by the Free Software Foundation.

# This source file implements the Coron_Roy_Vivek (CRV) / Pulkus_Vivek polynomial evaluation technique, with shared inputs. This code was written by Srinivas Vivek on 24-01-2016 13:50
# is modified here further otimize using SAT solvers.

# Evaluation of DES with 3 multiplications? Implementing and testing Jurgen Pulkus's method.
#A list containing the coefficients of polynomials for all the S-boxes are returned as integers. 
#The length of this list is 600 = 25*3*8.  A list containing evaluations of all the 72=3*3*8 linear polynomials 
#at every element of GF(2^8) is also generated. The irredducible polynomial used to represent  
#GF(2^8) is y^8 + y^4 + y^3 + y + 1. 
#    Last modified: 24-01-2016 13:50.

def polydecomp_sb_des(dom=6,bts=4,n=8,trm=2,cyls=[0,1,3,7]):
    # The domain is {0,1}^dom
    # Range is {0,1}^bts
    #Embedding field is 2^n
    # trm is the number of summands in the combining step
    # cyls is the list of pre-computed monomials input as a list of representatives of cyclotomic classes

    map_deg=[];
    for i in cyls:
        map_deg.append(i);
        if(i!= 2^n-1):
            j = i*2 % (2^n-1);
        else:
            j=i;
        while j!= i:
            map_deg.append(j);
            j = j*2 % (2^n-1);

# The code in the following lines to make 'map_deg'free of repetitions can be simplified using the 'set(map_deg)' command. 
    map_deg.sort();

    idl=len(map_deg);
    i=0;
    while (i < idl-1):
        if(map_deg[i]==map_deg[i+1]):
            map_deg.remove(map_deg[i]);
            idl -=1;
        else:
            i+=1;

    print ("L = ", map_deg)
    print ("|L|", idl)

    sb=[];

    #sb.append([1 % 2^bts for i in range(2^dom)]);     #Note: 2^dom instead of 2^n
    #sb.append([randint(0,2^bts-1) for i in range(2^dom)]);     #Note: 2^dom instead of 2^n

    sb.append([14, 0, 4, 15, 13, 7, 1, 4, 2, 14, 15, 2, 11, 13, 8, 1, 3, 10, 10, 6, 6, 12, 12, 11, 5, 9, 9, 5, 0, 3, 7, 8, 4, 15, 1, 12, 14, 8, 8, 2, 13, 4, 6, 9, 2, 1, 11, 7, 15, 5, 12, 11, 9, 3, 7, 14, 3, 10, 10, 0, 5, 6, 0, 13]);
    sb.append([15, 3, 1, 13, 8, 4, 14, 7, 6, 15, 11, 2, 3, 8, 4, 14, 9, 12, 7, 0, 2, 1, 13, 10, 12, 6, 0, 9, 5, 11, 10, 5, 0, 13, 14, 8, 7, 10, 11, 1, 10, 3, 4, 15, 13, 4, 1, 2, 5, 11, 8, 6, 12, 7, 6, 12, 9, 0, 3, 5, 2, 14, 15, 9]);
    sb.append([10, 13, 0, 7, 9, 0, 14, 9, 6, 3, 3, 4, 15, 6, 5, 10, 1, 2, 13, 8, 12, 5, 7, 14, 11, 12, 4, 11, 2, 15, 8, 1, 13, 1, 6, 10, 4, 13, 9, 0, 8, 6, 15, 9, 3, 8, 0, 7, 11, 4, 1, 15, 2, 14, 12, 3, 5, 11, 10, 5, 14, 2, 7, 12]);
    sb.append([7, 13, 13, 8, 14, 11, 3, 5, 0, 6, 6, 15, 9, 0, 10, 3, 1, 4, 2, 7, 8, 2, 5, 12, 11, 1, 12, 10, 4, 14, 15, 9, 10, 3, 6, 15, 9, 0, 0, 6, 12, 10, 11, 1, 7, 13, 13, 8, 15, 9, 1, 4, 3, 5, 14, 11, 5, 12, 2, 7, 8, 2, 4, 14]);
    sb.append([2, 14, 12, 11, 4, 2, 1, 12, 7, 4, 10, 7, 11, 13, 6, 1, 8, 5, 5, 0, 3, 15, 15, 10, 13, 3, 0, 9, 14, 8, 9, 6, 4, 11, 2, 8, 1, 12, 11, 7, 10, 1, 13, 14, 7, 2, 8, 13, 15, 6, 9, 15, 12, 0, 5, 9, 6, 10, 3, 4, 0, 5, 14, 3]);
    sb.append([12, 10, 1, 15, 10, 4, 15, 2, 9, 7, 2, 12, 6, 9, 8, 5, 0, 6, 13, 1, 3, 13, 4, 14, 14, 0, 7, 11, 5, 3, 11, 8, 9, 4, 14, 3, 15, 2, 5, 12, 2, 9, 8, 5, 12, 15, 3, 10, 7, 11, 0, 14, 4, 1, 10, 7,	1, 6, 13, 0, 11, 8, 6, 13]);
    sb.append([4, 13, 11, 0, 2, 11, 14, 7, 15, 4, 0, 9, 8, 1, 13, 10, 3, 14, 12, 3, 9, 5, 7, 12, 5, 2, 10, 15, 6, 8, 1, 6, 1, 6, 4, 11, 11, 13, 13, 8, 12, 1, 3, 4, 7, 10, 14, 7, 10, 9, 15, 5, 6, 0, 8, 15, 0, 14, 5, 2, 9, 3, 2, 12]);
    sb.append([13, 1, 2, 15, 8, 13, 4, 8, 6, 10, 15, 3, 11, 7, 1, 4, 10, 12, 9, 5, 3, 6, 14, 11, 5, 0, 0, 14, 12, 9, 7, 2, 7, 2, 11, 1, 4, 14, 1, 7, 9, 4, 12, 10, 14, 8, 2, 13, 0, 15, 6, 12, 10, 9, 13, 0, 15, 3, 3, 5, 5, 6, 8, 11]);
    
    F2 = GF(2)
    R.<y> = PolynomialRing(F2);

    irrd=[y^8+y^4+y^3+y+1];

    #for i in R.polynomials(of_degree = n):
    #if i.is_irreducible():
    #irrd.append(i);
    

    for i0 in [irrd[0]]:
        fincfsint=[]; # This list consists of all the coefficients of all the polynomials for all the S-boxes.
        linptab=[];   # This list consists of the polynomial evaluation table for all the (three) linear polynomials for each of the S-boxes.
        F.<a> = GF(2**n, name='a', modulus=i0);
        gfmap = [];
        for i1 in range(2^n): #n=8
            tmp1 = Integer(i1).bits();
            #print(tmp1)
            t1 = 0;
            t2 = 0;
            j1 = 0;
            for b7 in tmp1:
                t1 = t1 + b7 * a^j1;
                t2 = t2 + b7 * y^j1;
                j1 = j1+1;
            gfmap.append([t1,t2]);
        #print("gfmap=",gfmap);
        #print(F)
        x = PolynomialRing(F,(2*trm-1)*idl+1,'x',order='lex').gens(); # Looking for the decomposition sum(p_i*q_i)+p_trm
        w = PolynomialRing(F2,2*n+1,'w',order='lex').gens();
        #print("x=",x)

        def fnd_elm(tmp4):
            r1 = 0;
            while (gfmap[r1][0] != tmp4):
                r1 = r1 + 1;
            tmp6 = (gfmap[r1][1])*y^0;
            return [tmp6[i] for i in range(n)]

        def my_bits(z):     # Multiply two elements of F_{2^n} with unknown bit-cofficients, 
            #substitute latter one of them with the bit-vector z, and return a list of linear polynomials
            #for each bit of the product after modular reduction.
            p1 = sum([w[1+i]*(w[0]^i) for i in range(n)])*sum([w[n+1+i]*(w[0]^i) for i in range(n)]);
            tmp2 = [p1.coefficient({w[0]:i}) for i in range(2*n-1)]
            d2 = 2*(n-1)
            d3 = n
            for i3 in range(d2-d3+1):
                tmp3 = tmp2[d2-i3];
                for j3 in range(d3+1):
                    tmp2[d2-i3-j3] = tmp2[d2-i3-j3] + tmp3*i0[d3-j3];
            l0=[];
            for h2 in range(n):
                for h3 in range(n):
                    tmp2[h2] = tmp2[h2].subs({w[n+1+h3]:z[h3]})
                l0.append(tmp2[h2]);
            return l0;

        p0 = 0;
        for j in range(trm-1):
            p0 += sum([x[i+(2*j)*idl+1]*(x[0]^map_deg[i]) for i in range(idl)]) * sum([x[i+(2*j+1)*idl+1]*(x[0]^map_deg[i]) for i in range(idl)]);
        #Why python overflow must occur in the following command when n=16?
        
        p_sum = p0+ sum([x[i+2*(trm-1)*idl+1]*(x[0]^map_deg[i]) for i in range(idl)]);# P0*P1 + P2*P3 +...+ P_{2*trm-4}*P_{2*trm-3} +  P_{2*trm-2}		(P_{2*i+1} are choosen)
        #print(p_sum)
        RL.<u> = PolynomialRing(F,'u');
        A = matrix(F,2^dom * bts,trm*idl*n);    # Note: 2^dom * bts instead of 2^n * bts. #t.|l|.n
        cfsintb = [];
        plycfs=[];
        for j in range(trm-1):
            cfsintb.append( [randint(0,2^n-1) for i in range(idl)] );     # Coefficients, represented as integers, of the chosen polynomials P_{2*i+1}
            plycfs.append( [gfmap[i][0] for i in cfsintb[j]] );     # Coefficients, represented as field elements of F_{2^n}, of the chosen polynomials P_{2*i+1}
            
        
        tmp6 = p_sum;
        
         #---------------------------------------------------------------------------------
        #print(" *****************")
        for h2 in range(idl):#25
            for j in range(trm-1):#1
                x_end=(2*j+1*idl+1+h2)
        x_end=x_end+1
        
        print("x_end=",x_end)
        tmp_x=[];
        
        for r0 in range(2^dom):     #Note: 2^dom instead of 2^n
            tmp_prior = tmp6.subs({x[0]:gfmap[r0][0]});
            if tmp_prior.constant_coefficient() != 0:
                print ("Error")
            tmp_x.append(tmp_prior)
        #print("tmp_x=",tmp_x)
        l_end=(2*trm-1)*idl + 1
        print("l_end=",l_end)
        SatSolver = open('PATH TO YOUR SATSOLVER.EQS FILE', 'w')
        
        
        def replace_z(rep,st,z_assign):
            firstDelPos=st.find("(") 
            
            while(firstDelPos != -1):
                secondDelPos=st.find(")")
                #print(secondDelPos)
                rep_st="z"+str(rep)
                z_assign.append(rep)
                rep_sent=rep_st+"=1"
                #print("rep_sent ",rep_sent)
                rep=rep+1
                st=st.replace(st[firstDelPos:secondDelPos+1],rep_st)
                #print("st=",st)
                firstDelPos=st.find("(") 
            #print(st)    
            return(rep,st,z_assign)
        
        rep=1;
        z_assign=[]
        for xi in tmp_x:
            st=str(xi)
            #print("initially st=",st)
            #xi.subs((a+1):1)
            #tm="("+str(gfmap[])
            #st=st.replace("(a)","a")
            #st=st.replace("(a + 1)","z")
            
            #print("st=",st)
            rep,st,z_assign=replace_z(rep,st,z_assign)
            tmp_xst="tmp_x = "+str(st)
            SatSolver.write(str(tmp_xst))
            SatSolver.write(str("\n"))
            
        #print("finally st=",st)    
        #print("z_assign ",z_assign)
        for i in z_assign:
            z_st="z"+str(i)+"=1"
            SatSolver.write(str(z_st))
            SatSolver.write(str("\n"))
            
        SatSolver.write(str("tmp_x = 1"))
        SatSolver.write(str("\n"))
        #for j in range(1,trm+1):
         #   st2="x"
          #  st2=st2 + str(j) + " = 0"
           # SatSolver.write(str(st2))
            #SatSolver.write(str("\n"))
        for j in range(2*idl+1):
            if(j!=0):
                st2="x"
                st2=st2 + str(j) + " = 1"
                #print(st2)
                SatSolver.write(str(st2))
                SatSolver.write(str("\n"))
        
        SatSolver.close()
        inp="0"
        while(inp is "0"):
            inp=input("Press enter to indicate conversion is completed")
        print("conversion completed")
        
        file1 = open(r'PATH TO YOUR CLAIM FILE.txt','r')
        satConv=file1.read()
        #print(satConv)
        satConv=satConv.split("\n")
        #print(satConv)
        satOpt=[];
        for i in range(1,l_end):
            st="x"+str(i)+"=0"
            if(st in satConv):
                satOpt.append(i)
        #print("satopt=",satOpt)
        #print("tmp6=",tmp6)
        for i in satOpt:
            tmp6=tmp6.subs({x[i]:0});
        #print("tmp6=",tmp6)
        #print(" *****************")
        #----------------------------------------------------------------------------------
        
        for h2 in range(idl):#25
            for j in range(trm-1):#1
                tmp6 = tmp6.subs({x[(2*j+1)*idl+1+h2]:plycfs[j][h2]});
            
        
       
        for r0 in range(2^dom):     #Note: 2^dom instead of 2^n
            if r0%10 == 0:
                print ("r0 = ", r0)
            tmp7 = tmp6.subs({x[0]:gfmap[r0][0]});
            
            if tmp7.constant_coefficient() != 0:
                print ("Error")

            for c1 in range(idl):
                for c2 in range(trm):
                    tmp8 = my_bits(fnd_elm(tmp7.coefficient({x[(2*c2)*idl+1+c1]:1})))
                    for c3 in range(n):
                        for r1 in range(bts):
                            A[bts*r0+r1,n*(c2*idl+c1)+c3] = tmp8[r1].coefficient({w[1+c3]:1});

        print ("Rank",rank(A))
        #print list(A)

        def gf2int(z):
            tmp2=0;
            for i in range(len(z)):
                if(z[i] == 1):
                    tmp2 += 2^i
            return tmp2
    
        for tmpsb0 in sb:
            B = vector(F,2^dom*bts);     # Note: 2^dom * bts instead of 2^n * bts. Vector B has to be reset everytime.
#           fincfsint=[];         # This list consists of all the coefficients of all the polynomials for all the S-boxes.
            plys=[]; #Note: 2^dom instead of 2^n
            
        
            for r2 in range(2^dom):
                tmp5 = Integer(tmpsb0[r2]).bits();
                for h1 in range(len(tmp5)):
                    B[bts*r2+h1] = tmp5[h1];
    
            print ("Solving Sbox",sb.index(tmpsb0),", irrd index = ",irrd.index(i0),", irrd = ", i0);
            #print(B)
            try:
                X = A.solve_right(B);
            except(ValueError) as msg0:
                print(msg0);
                print("S-box",tmpsb0);
                print(plycfs);
                continue;

            slncfs = []
            cfsints = []
            for j in range(trm):
                slncfs.append([])
                cfsints.append([])
                for i in range(idl):
                    slncfs[j].append(gfmap[gf2int(X[n*(j*idl+i):n*((j*idl+i)+1)])][0])
                    cfsints[j].append(gf2int(X[n*(j*idl+i):n*((j*idl+i)+1)]))

            for j in range(trm-1):
                plys.append(sum([slncfs[j][i]*u^map_deg[i] for i in range(idl)] ));
                plys.append(sum([plycfs[j][i]*u^map_deg[i] for i in range(idl)] ));
            plys.append(sum([slncfs[trm-1][i]*u^map_deg[i] for i in range(idl)] ));	

            lp0 = 0
            for i in range(trm-1):
                lp0 += plys[2*i]*plys[2*i+1]
            lp0 += plys[2*trm-2];
            check=0
            for r3 in range(2^dom):     #Note: 2^dom instead of 2^n
                if(tmpsb0[r3] != ( Integer(gf2int(fnd_elm(lp0.substitute(u=gfmap[r3][0]))))%(2^bts)) ):
                    check=-1
                    print("Consistency Error",r3, tmpsb0[r3], ( Integer(gf2int(fnd_elm(lp0.substitute(u=gfmap[r3][0]))))%(2^bts)));
            if(check==0):
                print("same")
            #The following code is specific to DES S-boxes
            linplys=[];
            
            linplys.append(sum([slncfs[0][1]*u, slncfs[0][2]*u^2, slncfs[0][4]*u^4, slncfs[0][7]*u^8, slncfs[0][10]*u^16,  slncfs[0][13]*u^32, slncfs[0][16]*u^64,  slncfs[0][19]*u^128]));
            linplys.append(sum([slncfs[0][3]*u, slncfs[0][5]*u^2, slncfs[0][8]*u^4, slncfs[0][11]*u^8, slncfs[0][14]*u^16,  slncfs[0][17]*u^32, slncfs[0][22]*u^64,  slncfs[0][20]*u^128]));
            linplys.append(sum([slncfs[0][6]*u, slncfs[0][9]*u^2, slncfs[0][12]*u^4, slncfs[0][15]*u^8, slncfs[0][18]*u^16,  slncfs[0][24]*u^32, slncfs[0][23]*u^64,  slncfs[0][21]*u^128]));

            linplys.append(sum([plycfs[0][1]*u, plycfs[0][2]*u^2, plycfs[0][4]*u^4, plycfs[0][7]*u^8, plycfs[0][10]*u^16,  plycfs[0][13]*u^32, plycfs[0][16]*u^64,  plycfs[0][19]*u^128]));
            linplys.append(sum([plycfs[0][3]*u, plycfs[0][5]*u^2, plycfs[0][8]*u^4, plycfs[0][11]*u^8, plycfs[0][14]*u^16,  plycfs[0][17]*u^32, plycfs[0][22]*u^64,  plycfs[0][20]*u^128]));
            linplys.append(sum([plycfs[0][6]*u, plycfs[0][9]*u^2, plycfs[0][12]*u^4, plycfs[0][15]*u^8, plycfs[0][18]*u^16,  plycfs[0][24]*u^32, plycfs[0][23]*u^64,  plycfs[0][21]*u^128]));
            
            linplys.append(sum([slncfs[1][1]*u, slncfs[1][2]*u^2, slncfs[1][4]*u^4, slncfs[1][7]*u^8, slncfs[1][10]*u^16,  slncfs[1][13]*u^32, slncfs[1][16]*u^64,  slncfs[1][19]*u^128]));
            linplys.append(sum([slncfs[1][3]*u, slncfs[1][5]*u^2, slncfs[1][8]*u^4, slncfs[1][11]*u^8, slncfs[1][14]*u^16,  slncfs[1][17]*u^32, slncfs[1][22]*u^64,  slncfs[1][20]*u^128]));
            linplys.append(sum([slncfs[1][6]*u, slncfs[1][9]*u^2, slncfs[1][12]*u^4, slncfs[1][15]*u^8, slncfs[1][18]*u^16,  slncfs[1][24]*u^32, slncfs[1][23]*u^64,  slncfs[1][21]*u^128]));

            for i in linplys:
                for j in range(2^n):
                    linptab.append(Integer(gf2int(fnd_elm(i.substitute(u=gfmap[j][0])))));
            
            for i in cfsints[0]:
                fincfsint.append(i);
            for i in cfsintb[0]:
                fincfsint.append(i);
            for i in cfsints[1]:
                fincfsint.append(i);

            for r3 in range(2^dom):
                tmp9x=gfmap[r3][0];
                tmp9x3= Integer(gf2int(fnd_elm((tmp9x)^3)));
                tmp9x7= Integer(gf2int(fnd_elm((tmp9x)^7)));
            
            
            pp1= gfmap[ fincfsint[75*sb.index(tmpsb0)] ^^ linptab[2304* sb.index(tmpsb0)+r3] ^^ linptab[2304* sb.index(tmpsb0)+256+tmp9x3]  ^^ linptab[2304*sb.index(tmpsb0)+512+tmp9x7] ][0]
            pq1= gfmap[ fincfsint[75*sb.index(tmpsb0)+25] ^^ linptab[2304* sb.index(tmpsb0)+768+r3] ^^ linptab[2304* sb.index(tmpsb0)+768+256+tmp9x3]  ^^ linptab[2304*sb.index(tmpsb0)+768+512+tmp9x7] ][0]
            pp2= gfmap[ fincfsint[75*sb.index(tmpsb0)+25*2] ^^ linptab[2304* sb.index(tmpsb0)+768*2+r3] ^^ linptab[2304* sb.index(tmpsb0)+768*2+256+tmp9x3]  ^^ linptab[2304*sb.index(tmpsb0)+768*2+512+tmp9x7] ][0]
            
            if(lp0.substitute(u=gfmap[r3][0]) !=  (pp1*pq1+pp2)):
                print("Consistency Error",r3,lp0.substitute(u=gfmap[r3][0]),(pp1*pq1+pp2));
            


        target = open('des_tables.txt', 'w')
        target.write(str(fincfsint))
        #print(str(fincfsint))
        target.write("\n\n")
        target.write(str(linptab))

    return 0

polydecomp_sb_des()
