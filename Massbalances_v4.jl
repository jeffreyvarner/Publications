using Debug
using ODE
using Sundials
using Gadfly

#Alias the species in each compartment
#species[1]=concentration of A in heart
#species[2]=concentration of B in heart
#species[3]=concentration of A in aorta
#species[4]=concentration of B in aorta
#species[5]=concentration of A in organ 1
#species[6]=concentration of B in organ 1
#species[7]=concentration of A in organ 2
#species[8]=concentration of B in organ 2
#species[9]=concentration of A in vena cava
#species[10]=concentration of B in vena cava
#species[11]=concentration of A in wound space
#species[12]=concentration of B in wound space
#species[13]=Cardiac Output i.e. Q_heart_out
#species[14]=volume_blood_heart - volume[1]
#species[15]=volume_blood_aorta - volume[2]
#species[16]=volume_blood_organ1 - volume[3]
#species[17]=volume_blood_organ2 - volume[4]
#species[18]=volume_blood_venacava - volume[5]
#species[19]=volume_blood_wound - volume[6]


#-----------------------------------------------------------------
#Definitions of flow Rates when there is an injury to the aorta
#-----------------------------------------------------------------

#Flow rates in and out of heart
#Q_heart_out=1; #Flow rate out of the heart into aorta - aorta in the main artery. This is the cardiac output

#Assume there is a loss of blood due to injury to the aorta
#Q_out_aorta_wound1=eta[3]*Q_heart_out; #Flow rate out of blood vessel into the wound compartment at the time of injury
#Q_out_wound=eta[4]*Q_out_aorta_wound1; #Loss of blood from wound compartment
#Q_in_aorta_wound1=eta[5]*Q_out_aorta_wound1; #Flow of blood from wound compartment to aorta

#Flow into organ 1
#Q1_in = eta[1]*(Q_heart_out-Q_out_aorta_wound1+Q_in_aorta_wound1);
#Q1_out = eta[1]*(Q_heart_out-Q_out_aorta_wound1+Q_in_aorta_wound1);

#Flow into organ 2
#Q2_in = eta[2]*(Q_heart_out-Q_out_aorta_wound1+Q_in_aorta_wound1);
#Q2_out = eta[2]*(Q_heart_out-Q_out_aorta_wound1+Q_in_aorta_wound1);

#Q_heart_in=Q1_out+Q2_out;  #Flow rate into the heart from the inferior/superior vena cava - the main veins the carry deoxygenated blood into the heart

#---------------------------------------------------------------------------------

#Fraction of flow into each organ. Blood flow into each organ is going to be different. Assuming two organs for now. Last two eta's are fraction of flow back
#from wound compartment into the arterial pool and out of the wound compartment
#eta[1]- fraction of flow (cardiac output - blood loss) into organ 1
#eta[2]- fraction of flow (cardiac output - blood loss) into organ 2
#eta[3]- fraction of flow into wound
#eta[4]- fraction of wound flow outside - loss of blood
#eta[5]- fraction of wound flow back to aorta

Cardiac_output=1:10;
eta=[0.25,0.75,0.03,0.5,0.5];
#----------------------------------------------------
#Use organ volumes (dimensionless) rather than blood volumes for now
volume_heart=0.2;
volume_aorta=1;
volume_organ1=0.3;
volume_organ2=0.5;
volume_venacava=1;
volume_injury=0.007;

#----------------------------------------------------
initial_volumes=float([volume_heart,volume_aorta,volume_organ1,volume_organ2,volume_venacava,volume_injury]);
volume=float(initial_volumes);
DC=0.5; #dimensionless constant similar to Damkohler number
K_m=1.0; #Dimensionless Michaelis-Menten constant

function Massbalances(TSIM, species, speciesdot)

  #Doing -ve of outflux-influx  here
  #Massbalances in heart
  speciesdot[1]=-((species[13]*species[1])-((eta[1]+eta[2])*(1-eta[3]+eta[5]*eta[3])*species[13]*species[9]))/volume[1];
  speciesdot[2]=-((species[13]*species[2])-((eta[1]+eta[2])*(1-eta[3]+eta[5]*eta[3])*species[13]*species[10]))/volume[1];

  #Massbalances in aorta
  speciesdot[3]=-((eta[1]+eta[2])*(1-eta[3]+eta[5]*eta[3])*species[13]*species[3]+(eta[3]*species[13]*species[3])-(species[13]*species[1]+eta[5]*eta[3]*species[13]*species[11]))/volume[2];
  speciesdot[4]=-((eta[1]+eta[2])*(1-eta[3]+eta[5]*eta[3])*species[13]*species[4]+(eta[3]*species[13]*species[4])-(species[13]*species[2]+eta[5]*eta[3]*species[13]*species[12]))/volume[2];

  #Massbalances in organ1
  speciesdot[5]=-(eta[1]*(1-eta[3]+eta[5]*eta[3])*species[13])*(species[5]-species[3])/volume[3];
  speciesdot[6]=-(eta[1]*(1-eta[3]+eta[5]*eta[3])*species[13])*(species[6]-species[4])/volume[3];

  #Massbalances in organ2
  speciesdot[7]=-(eta[2]*(1-eta[3]+eta[5]*eta[3])*species[13])*(species[7]-species[3])/volume[4];
  speciesdot[8]=-(eta[2]*(1-eta[3]+eta[5]*eta[3])*species[13])*(species[8]-species[4])/volume[4];

  #Massbalances in venacava
  speciesdot[9]=-((((eta[1]+eta[2])*(1-eta[3]+eta[5]*eta[3])*species[13])*species[9])-((eta[1]*(1-eta[3]+eta[5]*eta[3])*species[13])*species[5]+((eta[2]*(1-eta[3]+eta[5]*eta[3])*species[13])*species[7])))/volume[5];
  speciesdot[10]=-((((eta[1]+eta[2])*(1-eta[3]+eta[5]*eta[3])*species[13])*species[10])-((eta[1]*(1-eta[3]+eta[5]*eta[3])*species[13])*species[6]+((eta[2]*(1-eta[3]+eta[5]*eta[3])*species[13])*species[8])))/volume[5];

  #Massbalances in wound compartment
  speciesdot[11]=-(((eta[4]*eta[3]+eta[5]*eta[3])*species[13]*species[11]-(eta[3]*species[13])*species[3])/volume[6])-(DC)*(species[11])/(K_m+species[11]);
  speciesdot[12]=-(((eta[4]*eta[3]+eta[5]*eta[3])*species[13]*species[12]-(eta[3]*species[13])*species[4])/volume[6])+(DC)*(species[11])/(K_m+species[11]);

  #Flow rates as diffeqs
  speciesdot[13]= -eta[4]*eta[3]*species[13]; #Flow out of wound will be proportional to caridiac output
end

#------------------SOLVE THE MODEL WITH INJURY THAT CAUSES BLOOD LOSS-----------------------------
#z=zeros(,100);
for i in Cardiac_output
  IC=float([1,0,1,0,1,0,1,0,1,0,0,0,Cardiac_output[i]]);
  TSIM = float([1:1:100]);
  XSIM = Sundials.cvode(Massbalances, IC, TSIM);
  x=TSIM;
  y=XSIM[:,12];
  s=string("./results",i,".txt");
  writedlm(s,y);
end
