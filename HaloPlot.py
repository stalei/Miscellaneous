import matplotlib.pyplot as plt
import numpy as np
from matplotlib.legend_handler import HandlerLine2D
#import seaborn as sns
import tkinter

########################################################### Data source
# #id num_p mvir mbound_vir rvir vmax rvmax vrms x y z vx vy vz Jx Jy Jz E Spin PosUncertainty VelUncertainty
#bulk_vx bulk_vy bulk_vz BulkVelUnc n_core m200b m200c m500c m2500c Xoff Voff spin_bullock b_to_a c_to_a A[x] A[y] A[z]
#b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c) A[z](500c) Rs Rs_Klypin T/|U| M_pe_Behroozi M_pe_Diemer 
#Halfmass_Radius idx i_so i_ph num_cp mmetric








def PlotHalos(event):
	data=np.genfromtxt('halos_0.0.ascii', skip_header=18)#,names=True, skip_header=5)

	##Read parameters

	id=np.array(data[:,0])
	count= len(id)

	num_p=np.array(data[:,1])
	mvir=np.array(data[:,2])
	rvir=np.array(data[:,4])
	x=np.array(data[:,8])
	y=np.array(data[:,9])
	z=np.array(data[:,10])


	## Criterias

	MHL=1.0e11
	DHL=1.0e8
	ProbeRadius=0.3 # kpc -> Mpc
	NumLimit=40
	nbins=30
	nbins2=50
	## Mass Sample

	#Mglog=MglogAll[NumG>NumLimit]



	#MMassiveG=Mg[Mg>MHL]


	print("Min # of particles in a halo is:", NumLimit)

	print("No of halos above", "%.4g"%MHL,"is:")
	#print "Gadget:",len(MMassiveG)


	#MDwarfG=MG2[NumDwarfG>NumLimit]

	print("We selected this massive halo:")
	#print "Gadget=%.4g"%MMassiveG[MMG_index],"solar mass & coordinates:",XMassiveG[MMG_index],",",YMassiveG[MMG_index],",",ZMassiveG[MMG_index]

	#print "No of sattelites closer than",ProbeRadius*1000," kpc is:"
	#print "Gadget:",len(rSampleG)


'''
########################################################## plots
fig = plt.figure(1)
fig.suptitle('CoSANG vs N-Body ')

ax1 = fig.add_subplot(221)
ax1.set_xlabel('$Log (M_{halo})$')
ax1.set_ylabel('$N$')
ax1.set_title('Dwarf Halos Mass aboundance')
#ax1.plot(d1['star_age'],d1['center h1']) #plot of main data
ax1.hist(np.log10(MDwarfG),bins=nbins, log=False, histtype='step', alpha=0.9,color='blue',label='Gadget')
ax1.hist(np.log10(MDwarfC0),bins=nbins,log=False, histtype='step', alpha=0.9,color='black',label='NoDisk')
ax1.hist(np.log10(MDwarfC),bins=nbins,log=False, histtype='step', alpha=0.9,color='green',label='CoSANG')

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



ng, binsg, patchesg= ax3.hist(np.log10(rSampleG),bins=nbins,log=False,cumulative=True, histtype='step', alpha=0.8,color='blue',label='Gadget')
nc0, binsc0, patchesc0= ax3.hist(np.log10(rSampleC0),bins=nbins,log=False,cumulative=True,histtype='step', alpha=0.8,color='black',label='NoDisk')
nc, binsc, patchesc= ax3.hist(np.log10(rSampleC),bins=nbins,log=False, cumulative=True,histtype='step', alpha=0.8,color='green',label='CoSANG')
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
ng2=ng[ng>0]
nc2=nc[ng>0]
nc02=nc0[ng>0]
line1,= ax4.plot(nc2/ng2,'k+',label='nc/ng')
line2,= ax4.plot(nc02/ng2,'ko', label='nc0/ng')
#ax4.plot(np.log10(ng2/nc2),'r')

ax4.legend(handler_map={line1: HandlerLine2D(numpoints=4)})

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
'''






##Start GUI
window = tkinter.Tk()
# to rename the title of the window
window.title("Halo Plot")
# pack is used to show the object in the window
label_top = tkinter.Label(window, text = "Halo Plot is desined to overview halo information using Rockstar output",fg="white", bg="black").pack(fill = "x")


# creating 2 frames TOP and BOTTOM
left_frame = tkinter.Frame(window).pack()
right_frame = tkinter.Frame(window).pack(side = "right")

label_left = tkinter.Label(left_frame, text = "Plot").pack()



label_right = tkinter.Label(right_frame, text = "Controls").pack()


btn = tkinter.Button(right_frame, text = "Plot")
btn.bind("<Button-1>", PlotHalos) # 'bind' takes 2 parameters 1st is 'event' 2nd is 'function'
btn.pack()



label_count = tkinter.Label(left_frame, textvariable =count).pack()


window.mainloop()