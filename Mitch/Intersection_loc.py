
from shapely.geometry import LineString
from EpEq_space import * 

'''
Calculates the intersection point for 2-sigma containment for electron recoils and nuclear recoils in EQ_EP space. 
'''





def Intersection_loc(DF_er,DF_nr,cut_idx,bins,fano):


    EP_er_Center, EQ_er_edge, EP_nr_Center, EQ_nr_edge = EpEQ_bin(DF_er,DF_nr,cut_idx,bins) #binned 

    

    m1,c1,r1,p1,se1= stats.linregress(EP_er_Center,EQ_er_edge)
    m2,c2,r2,p2,se2= stats.linregress(EP_nr_Center,EQ_nr_edge)


    EQ_er_new = []
    for x in EP_er_Center: 
        new1 = x*m1 +c1
        EQ_er_new.append(new1) 

    EQ_nr_new = []
    for x in EP_nr_Center: 
        new2 = x*m2 +c2
        EQ_nr_new.append(new2) 

    '''Adds y-intercept to data array'''
    EP_er_Center.append(0)
    EQ_er_edge.append(c1)
    EQ_er_new.append(c1)

    EP_nr_Center.append(0)
    EQ_nr_edge.append(c2)
    EQ_nr_new.append(c2)
    


    w = np.sort(EP_nr_Center)
    x = np.sort(EQ_nr_new)
    y = np.sort(EP_er_Center) 
    z = np.sort(EQ_er_new)
    idx = len(w)-1
    '''Finds Intersection of fitted lines'''
    line1 = LineString([(w[0],x[0]),(w[idx],x[idx])])
    line2 = LineString([(y[0],z[0]),(y[idx],z[idx])])

    Intersection = line2.intersection(line1)




    '''This plot showes the Ep_Eq space for both electron and nuclear recoils. '''
    plt.figure(figsize=(9.0,8.0))

    '''Data'''
    plt.scatter(DF_er.EP,DF_er.EQ,color = 'black',s= 2.5,label = 'Electron Recoils')
    plt.scatter(DF_nr.EP,DF_nr.EQ,s= 2.5,label = 'Nuclear Recoils')
    '''Fitted Points From Distribution'''
    plt.scatter(EP_er_Center,EQ_er_edge,color = 'red') 
    plt.scatter(EP_nr_Center,EQ_nr_edge,color = 'black')
    '''Fitted line'''   
    plt.plot(EP_er_Center, EQ_er_new,'k--')
    plt.plot(EP_nr_Center, EQ_nr_new,'b--')


    plt.title('$E_P$ vs. $E_Q$ Fano = ' +str(fano), size = 16)
    plt.xlabel('$E_P$ [keV]',size = 14)
    plt.ylabel('$E_Q$ [keV]',size = 14)

    plt.grid(linestyle='-', linewidth=0.35)
    plt.legend()
    plt.savefig('/Users/Mitch 1/Desktop/Thesis_Plots/Cross_Over/Fitted_Fano = ' +str(fano)+'.png')

    plt.show()



    ''' This plot shows location of intersection'''


    plt.figure(figsize=(9.0,8.0))
    plt.plot(EP_er_Center, EQ_er_new,'k--')
    plt.plot(EP_nr_Center, EQ_nr_new,'b--')
    

    plt.scatter(DF_er.EP,DF_er.EQ,color = 'black',s= 2.5,label = 'Electron Recoils')
    plt.scatter(DF_nr.EP,DF_nr.EQ,s = 2.5,label = 'Nuclear Recoils')
    plt.plot(Intersection.x,Intersection.y,'rX',markersize = 12,label = 'Intersection')

    plt.scatter(EP_er_Center,EQ_er_edge) 
    plt.scatter(EP_nr_Center,EQ_nr_edge,color = 'black')

    plt.title('$E_P$ vs. $E_Q$ Fano = ' +str(fano) ,size = 16)
    plt.xlabel('$E_P$ [keV]',size = 14)
    plt.ylabel('$E_Q$ [keV]',size = 14)


    plt.xlim(-1,20)
    plt.ylim(-1,10)
    plt.grid(linestyle='-', linewidth=0.35)
    plt.legend(loc = 2)
    plt.savefig('/Users/Mitch 1/Desktop/Thesis_Plots/Cross_Over/Intersection_Fano = ' +str(fano)+'.png')

    plt.show()
    

    


    return Intersection 