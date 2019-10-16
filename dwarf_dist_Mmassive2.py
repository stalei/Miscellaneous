import matplotlib.pyplot as plt
import numpy as np
from matplotlib.legend_handler import HandlerLine2D
import seaborn as sns

########################################################### Data source
# logy for subhalo distplot
# MWish MMasive

dg=np.genfromtxt('halos_0.0_G.ascii', skip_header=18)#,names=True, skip_header=5)
dc=np.genfromtxt('halos_0.0_C.ascii',skip_header=18)#names=True, skip_header=5)
dc0=np.genfromtxt('halos_0.0_C0.ascii',skip_header=18)#names=True, skip_header=5)


IDG=np.array(dg[:,0])
IDC=np.array(dc[:,0])
IDC0=np.array(dc0[:,0])

NumG=np.array(dg[:,1])
NumC=np.array(dc[:,1])
NumC0=np.array(dc0[:,1])


#print len(dc[:,2])-len(dg[:,2]) # mvir
Mg=np.array(dg[:,2])
Mc=np.array(dc[:,2])
Mc0=np.array(dc0[:,2])

MglogAll=np.array(np.log10(dg[:,2]))
MclogAll=np.array(np.log10(dc[:,2]))
Mc0logAll=np.array(np.log10(dc0[:,2]))

#MglogAll=np.array(dg[:,2])
#MclogAll=np.array(dc[:,2])
#Mc0logAll=np.array(dc0[:,2])


Xg=np.array(dg[:,8])
Yg=np.array(dg[:,9])
Zg=np.array(dg[:,10])

Xc=np.array(dc[:,8])
Yc=np.array(dc[:,9])
Zc=np.array(dc[:,10])

Xc0=np.array(dc0[:,8])
Yc0=np.array(dc0[:,9])
Zc0=np.array(dc0[:,10])


RVirg=np.array(dg[:,4])
RVirc=np.array(dc[:,4])
RVirc0=np.array(dc0[:,4])

## Criterias

MHL=1.0e12
DHL=1.0e8
ProbeRadius=0.3 # kpc -> Mpc
NumLimit=100
nbins=30
nbins2=50
## Mass Sample

Mglog=MglogAll[NumG>NumLimit]
Mclog=MclogAll[NumC>NumLimit]
Mc0log=Mc0logAll[NumC0>NumLimit]


MMassiveG=Mg[Mg>MHL]
XMassiveG=Xg[Mg>MHL]
YMassiveG=Yg[Mg>MHL]
ZMassiveG=Zg[Mg>MHL]
RvirMassiveG=RVirg[Mg>MHL]

MMassiveC=Mc[Mc>MHL]
XMassiveC=Xc[Mc>MHL]
YMassiveC=Yc[Mc>MHL]
ZMassiveC=Zc[Mc>MHL]
RvirMassiveC=RVirc[Mc>MHL]

MMassiveC0=Mc0[Mc0>MHL]
XMassiveC0=Xc0[Mc0>MHL]
YMassiveC0=Yc0[Mc0>MHL]
ZMassiveC0=Zc0[Mc0>MHL]
RvirMassiveC0=RVirc0[Mc0>MHL]



rMassiveG=np.sqrt(XMassiveG**2.+YMassiveG**2.+ZMassiveG**2.)
rMassiveC=np.sqrt(XMassiveC**2.+YMassiveC**2.+ZMassiveC**2.)
rMassiveC0=np.sqrt(XMassiveC0**2.+YMassiveC0**2.+ZMassiveC0**2.)

print("Min # of particles in a halo is:", NumLimit )

print("No of halos above", "%.4g"%MHL,"is:")
print("Gadget:",len(MMassiveG))
print("NoDisk:",len(MMassiveC0))
print("CoSANG:",len(MMassiveC))


#print XMassiveG

NumDwarfG=NumG[Mg<DHL]
NumDwarfC0=NumC0[Mc0<DHL]
NumDwarfC=NumC[Mc<DHL]

Mg2=Mg[Mg<DHL]
print(len(Mg2))
MDwarfG_test=Mg2[NumDwarfG>NumLimit]
print(len(MDwarfG_test))
###############

MG2=Mg[Mg<DHL]
XG2=Xg[Mg<DHL]
YG2=Yg[Mg<DHL]
ZG2=Zg[Mg<DHL]
RVirG2=RVirg[Mg<DHL]

MC2=Mc[Mc<DHL]	
XC2=Xc[Mc<DHL]
YC2=Yc[Mc<DHL]
ZC2=Zc[Mc<DHL]
RVirC2=RVirc[Mc<DHL]


MC02=Mc0[Mc0<DHL]
XC02=Xc0[Mc0<DHL]
YC02=Yc0[Mc0<DHL]
ZC02=Zc0[Mc0<DHL]
RVirC02=RVirc0[Mc0<DHL]


###############


MDwarfG=MG2[NumDwarfG>NumLimit]
XDwarfG=XG2[NumDwarfG>NumLimit]
YDwarfG=YG2[NumDwarfG>NumLimit]
ZDwarfG=ZG2[NumDwarfG>NumLimit]
RvirDwarfG=RVirG2[NumDwarfG>NumLimit]

MDwarfC=MC2[NumDwarfC>NumLimit]
XDwarfC=XC2[NumDwarfC>NumLimit]
YDwarfC=YC2[NumDwarfC>NumLimit]
ZDwarfC=ZC2[NumDwarfC>NumLimit]
RvirDwarfC=RVirC2[NumDwarfC>NumLimit]



MDwarfC0=MC02[NumDwarfC0>NumLimit]
XDwarfC0=XC02[NumDwarfC0>NumLimit]
YDwarfC0=YC02[NumDwarfC0>NumLimit]
ZDwarfC0=ZC02[NumDwarfC0>NumLimit]
RvirDwarfC0=RVirC02[NumDwarfC0>NumLimit]


rDwarfG=np.sqrt(XDwarfG**2.+YDwarfG**2.+ZDwarfG**2.)
rDwarfC=np.sqrt(XDwarfC**2.+YDwarfC**2.+ZDwarfC**2.)
rDwarfC0=np.sqrt(XDwarfC0**2.+YDwarfC0**2.+ZDwarfC0**2.)


#print len(MDwarfG)
print("No of halos belove %.4g"%DHL,"is:")
print("Gadget:",len(MDwarfG))
print("NoDisk:",len(MDwarfC0))
print("CoSANG:",len(MDwarfC))



#### Find most massive halo
MMG_index=np.argmin(MMassiveG)
MMC_index=np.argmin(MMassiveC)
MMC0_index=np.argmin(MMassiveC0)

rMg=np.sqrt(XMassiveG[MMG_index]**2. + YMassiveG[MMG_index]**2. + ZMassiveG[MMG_index]**2.)
rMc0=np.sqrt(XMassiveC0[MMC0_index]**2. + YMassiveC0[MMC0_index]**2. + ZMassiveC0[MMC0_index]**2.)
rMc=np.sqrt(XMassiveC[MMC_index]**2. + YMassiveC[MMC_index]**2. + ZMassiveC[MMC_index]**2.)


print("We selected this massive halo:")
print("Gadget=%.4g"%MMassiveG[MMG_index],"solar mass & coordinates:",XMassiveG[MMG_index],",",YMassiveG[MMG_index],",",ZMassiveG[MMG_index])
print("NoDisk=%.4g"%MMassiveC0[MMC0_index],"solar mass & coordinates:",XMassiveC0[MMC0_index],",",YMassiveC0[MMC0_index],",",ZMassiveC0[MMC0_index])
print("CoSANG=%.4g"%MMassiveC[MMC_index],"solar mass & coordinates:",XMassiveC[MMC_index],",",YMassiveC[MMC_index],",",ZMassiveC[MMC_index])


delrG=np.sqrt((XDwarfG-XMassiveG[MMG_index])**2.+(YDwarfG-YMassiveG[MMG_index])**2.+(ZDwarfG-ZMassiveG[MMG_index])**2.)
delrC=np.sqrt((XDwarfC-XMassiveC[MMC_index])**2.+(YDwarfC-YMassiveC[MMC_index])**2.+(ZDwarfC-ZMassiveC[MMC_index])**2.)
delrC0=np.sqrt((XDwarfC0-XMassiveC0[MMC0_index])**2.+(YDwarfC0-YMassiveC0[MMC0_index])**2.+(ZDwarfC0-ZMassiveC0[MMC0_index])**2.)

rSampleG=1000.0*delrG[delrG<ProbeRadius]
rSampleC=1000.0*delrC[delrC<ProbeRadius]
rSampleC0=1000.0*delrC0[delrC0<ProbeRadius]

#massSampleG=MDwarfG[delrG<ProbeRadius]
#massSampleC=MDwarfG[delrC<ProbeRadius]
#massSampleC0=MDwarfG[delrC0<ProbeRadius]

print("Sample size:")
#print(len(rSampleG),len(massSampleG))
print(len(rSampleG),len(delrG))



#rSampleG=[]
#rSampleC=[]
#rSampleC0=[]

#for rmg in rDwarfG:
#	delrG=np.abs(rMg-rmg)
#	if delrG<ProbeRadius:
#		rSampleG.append(1000.*delrG)


#for rmc in rDwarfC:
#	delrC=np.abs(rMc-rmc)
#	if delrC<ProbeRadius:
#		rSampleC.append(1000.*delrC)

#delrG=np.abs(rMassiveG-rDwarfG)
#delrC=np.abs(rMassiveC-rDwarfC)

#rSampleG=delr[delrG<ProbeRadius]
#rSampleC=delr[delrC<ProbeRadius]

print("No of sattelites closer than",ProbeRadius*1000," kpc is:")
print("Gadget:",len(rSampleG))
print("NoDisk:",len(rSampleC0))
print("CoSANG:",len(rSampleC))

#print len(rSampleC)


########################################################## plots
fig = plt.figure(1)
fig.suptitle('CoSANG vs N-Body ')

ax1 = fig.add_subplot(221)
ax1.set_xlabel('$Log (M_{halo})$')
ax1.set_ylabel('$N$')
ax1.set_title('Dwarf Halos Mass aboundance')
#ax1.plot(d1['star_age'],d1['center h1']) #plot of main data
ax1.hist(np.log10(MDwarfG),linewidth=2, bins=nbins, log=False, histtype='step', alpha=0.9,color='blue',label='Gadget')
ax1.hist(np.log10(MDwarfC0),linewidth=2,bins=nbins,log=False, histtype='step', alpha=0.9,color='black',label='NoDisk')
ax1.hist(np.log10(MDwarfC),linewidth=2,bins=nbins,log=False, histtype='step', alpha=0.9,color='green',label='CoSANG')

ax1.legend(loc=2)
#line1, = ax1.plot(np.log10(d1['rvir']),d1['mvir'],linestyle='--', label='7$M_{\odot}$') # first extra plot for legend

#line2, = ax1.plot(np.log10(d2['star_age']),d2['log_R'],linestyle=':',label='8$M_{\odot}$') # second extra plot for legend

#line3, = ax1.plot(np.log10(d3['star_age']),d3['log_R'],linestyle='-.',label='9$M_{\odot}$') # second extra plot for legend
#line4, = ax1.plot(np.log10(d4['star_age']),d4['log_R'],linestyle='-',label='11$M_{\odot}$') 

#line5, = ax1.plot(np.log10(d1['star_age']),d1['log_cntr_T'],linestyle='--', label='7$M_{\odot}$')


#ax1.legend(handler_map={line1: HandlerLine2D(numpoints=4)})
 # first extra plot for legend

######################################################### 
ax2 = fig.add_subplot(222)

ax2.set_xlabel('$Log (M_{halo})$')
ax2.set_ylabel('$N$')
ax2.set_title('Mass aboundance')
colors=['green','black','blue']
labels=['CoSANG','NoDisk','Gadget']

ax2.hist([Mclog, Mc0log, Mglog], bins=nbins, log=True, color=colors,label=labels)
ax2.legend()
#########################################################
ax3 = fig.add_subplot(223)
ax3.set_xlabel('$log(d_{kpc})$')
ax3.set_ylabel('$N$')
ax3.set_title('Subhalo Distribution for ~%.1g halo'%MMassiveG[MMG_index])

#Mc=np.log10(dc[:,2])
#Mg=np.log10(dg[:,2])



ng, binsg, patchesg= ax3.hist(np.log10(rSampleG),linewidth=2,bins=nbins,log=False,cumulative=True, histtype='step', alpha=0.8,color='blue',label='Gadget')
nc0, binsc0, patchesc0= ax3.hist(np.log10(rSampleC0),linewidth=2,bins=nbins,log=False,cumulative=True,histtype='step', alpha=0.8,color='black',label='NoDisk')
nc, binsc, patchesc= ax3.hist(np.log10(rSampleC),linewidth=2,bins=nbins,log=False, cumulative=True,histtype='step', alpha=0.8,color='green',label='CoSANG')
#ax3.legend()


#ng, binsg, patchesg= ax3.hist(rSampleG,bins=100,log=True, alpha=0.6,color='green',label='NoDisk')
#nc, binsc, patchesc= ax3.hist(rSampleC,bins=100,log=True, alpha=0.4,color='blue',label='CoSANG')
#ax3.legend()
ax3.legend( loc =2)


#sns.distplot(np.log10(rSampleC), hist = False, kde = True,kde_kws = {'linewidth': 3},label = 'CoSANG')
#sns.distplot(np.log10(rSampleG), hist = False, kde = True,kde_kws = {'linewidth': 3},label = 'NoDisk')

#######################################################

ax4=fig.add_subplot(224)
ax4.set_xlabel('$Bin number~d_{kpc}$')
ax4.set_ylabel('$n_i/n_g$')
ax4.set_title('Subhalo Distribution')
colors=['grey','black']
labels=['CoSANG','NoDisk']

#binsc[binsc%2==1]
#ax4_x=np.linspace(0,ProbeRadius*1000,400)
#ax4.plot( ax4_x, nc/ng, '+')
#ax4.plot( nc/ng, 'k.')
#ng2=ng[ng>0]
#nc2=nc[ng>0]
#nc02=nc0[ng>0]
#line1,= ax4.plot(nc2/ng2,'k+',label='nc/ng')
#line2,= ax4.plot(nc02/ng2,'ko', label='nc0/ng')
#ax4.plot(np.log10(ng2/nc2),'r')

#line1,=ax4.plot(np.log10(rSampleG),np.log10(massSampleG) ,linewidth=2)

#ax4.legend(handler_map={line1: HandlerLine2D(numpoints=4)})

#######################################################
fig2 = plt.figure(2)
#fig2.set_title('Mass distribution')
#sns.kdeplot(Mc)
ax2_1=fig2.add_subplot(111)
ax2_1.set_xlabel('$Log(M)$')
ax2_1.set_ylabel('$n$')
ax2_1.set_title('Dwarf Halos Mass Distribution')
sns.distplot(np.log10(MDwarfC), hist = False, kde = True,kde_kws = {'linewidth': 3},label = 'CoSANG')
sns.distplot(np.log10(MDwarfC0), hist = False, kde = True,kde_kws = {'linewidth': 3},label = 'NoDisk')
sns.distplot(np.log10(MDwarfG), hist = False, kde = True,kde_kws = {'linewidth': 3},label = 'Gadget')
#fig2.legend(prop={'size': 16}, title = 'Airline')

###########################
plt.show()
