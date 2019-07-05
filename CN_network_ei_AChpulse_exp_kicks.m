function [ConMat,spktime,spkcell] = CN_network_ei_AChpulse_exp_kicks(wee,wii,wei,wie,rate,amp)
%[ConMat,spktime,spkcell] = CN_network_ei_AChpulse_exp_kicks(0.004,0.016,0.002,0.003,11,10);

N=1000; %total number of cells
num_ex=N*0.8; %number of excitatory cells
num_in=N*0.2; %number of inhibitory cells

ConMat=ConEI(num_ex,num_in,0.05,0.3,0.3,0.3); %generate connectivity matrix
W=[wee*ones(num_ex),wei*ones(num_ex,num_in);...
    wie*ones(num_in,num_ex),wii*ones(num_in)];     
ConMatW=ConMat.*W; %connectivity matrix with weight

t_all=4000; %total simulation time
t_pulse=2000; %time when ACh is applied
dt=0.05; %time step for numerical method
totalpts=t_all/dt;

gkse=0.6; %initial slow potassium channel conductance for E cells
gksi=0; %initial slow potassium channel conductance for I cells
t_increase=100; %ACh pulse rising time
t_decrease=1200; %ACh pulse decaying time

Iapp_ex=2.8140+0.6129*rand(1,num_ex); % baseline current for E cells
Iapp_in=-0.22+0.04*rand(1,num_in); % baseline current for I cells
Iapp=[Iapp_ex,Iapp_in];
Iapp_ini=[zeros(1,num_ex),0.57448+0.2607*rand(1,num_in)];
%baseline current for initialization

%setting parameters for neuron model
c=1;
gna=24;
gkdr=3;
gl=0.02;
ena=55;
ek=-90;
el=-60;

%setting parameters for synapses model
esyn_ex=0;
esyn_in=-75;
tau_d_ex=3;
tau_d_in=5.5;
tau_r=0.2;

%functions for neuron model
function h_inf=hinf(x)
    h_inf=1./(1+exp((x+53)./7));
end        

function tau_h=tauh(x)
	tau_h=0.37+2.78./(1+exp((x+40.5)./6));
end

function dh=dhdt(v,h)
    dh=(hinf(v)-h)./(tauh(v));
end

function n_inf=ninf(x)
	n_inf=1./(1+exp((-x-30)./10));
end
function tau_n=taun(x)
	tau_n=0.37+1.85./(1+exp((x+27)./15));
end
function dn=dndt(v,n)        
    dn=(ninf(v)-n)./(taun(v));
end

function z_inf=zinf(x)
	z_inf=1./(1+exp((-x-39)./5));
end

function dz=dzdt(v,z)
    tauz=75;
    dz=(zinf(v)-z)./tauz;
end

function m_inf=minf(x)
	m_inf=1./(1+exp((-x-30)./9.5));
end

%function for synaptic current
function I_syn=Isyn(v,t,spk)
    Ioutput_single_term1=zeros(1,N);
    Ioutput_single_term2=zeros(1,N);
    for k=1:num_ex
        spikepre=spk(:,k);
        spikepre(spikepre==0)=[];
        if isempty(spikepre)==0 
            Ioutput_single_term1(k)=sum(exp(-(t-spikepre)/tau_d_ex)-...
                exp(-(t-spikepre)/tau_r));
            Ioutput_single_term2(k)=sum((exp(-(t-spikepre)/tau_d_ex)-...
                exp(-(t-spikepre)/tau_r))*esyn_ex);
        end
    end
    for k=(num_ex+1):N
        spikepre=spk(:,k);
        spikepre(spikepre==0)=[];
        if isempty(spikepre)==0 
            Ioutput_single_term1(k)=sum(exp(-(t-spikepre)/tau_d_in)-...
                exp(-(t-spikepre)/tau_r));
            Ioutput_single_term2(k)=sum((exp(-(t-spikepre)/tau_d_in)-...
                exp(-(t-spikepre)/tau_r))*esyn_in);
        end
    end
    I_syn=(Ioutput_single_term1*ConMatW).*v-Ioutput_single_term2*ConMatW;
end

%function for ACh pulse
function gks=gks_pulse(t)
    gkse_t=gkse;
    if t>=t_pulse && t<t_pulse+t_increase 
        gkse_t=gkse-gkse*(t-t_pulse)/t_increase;
    else
        if t>=t_pulse+t_increase
            gkse_t=gkse-gkse*exp(-(t-t_pulse-t_increase)/(t_decrease/3));
        end          
    end
    gks=[gkse_t*ones(1,num_ex),gksi*ones(1,num_in)];
end

%generate random Poisson noise
Kick_time=[];
Kick_cell=[];
const=300;
for index=1:N
    temp_exp=exprnd(1/rate,const,1)*1000;
    Kick_time_temp=(tril(ones(const))*temp_exp)';
    Kick_time_temp(Kick_time_temp>t_all)=[];
    Kick_cell_temp=index*ones(1,length(Kick_time_temp));
    Kick_time=[Kick_time,Kick_time_temp];
    Kick_cell=[Kick_cell,Kick_cell_temp];
end

%main function for neuron model
function dv=dvdt(v,h,n,z,t,spk)
    if t>100    
        I_syn=Isyn(v,t,spk);
    else
        I_syn=0;
    end
    if t<=100
        Iapp_add=Iapp_ini;
    else
         Iapp_add=0;
    end
    K=zeros(1,N);
    id=find((t-Kick_time)>=0 & (t-Kick_time)<=1);
    cell=Kick_cell(id);
    K(cell)=1;
    gks=gks_pulse(t);
    Iapp_kick=amp*K.*ones(1,N);
    dv=(Iapp+Iapp_add+Iapp_kick-I_syn-gna*(minf(v).^3).*h.*(v-ena)-...
        gkdr*(n.^4).*(v-ek)-gks.*z.*(v-ek)-gl*(v-el))/c;    
end

%Initialize membrane voltages and gating variables
v=-62+40*rand(1,N);
h=0.2+0.6*rand(1,N);
n=0.2+0.6*rand(1,N);
z=0.15+0.1*rand(1,N);
ind=1;

spktime_size_pre=(t_all*N*250)/1000;
spktime=zeros(1,spktime_size_pre);
spkcell=zeros(1,spktime_size_pre);
spk=zeros(3,N);

vt2=v;

%solve by RK4
for i=1:totalpts
    t=dt*(i-1);
    vt1=vt2;
    vt2=v;
    vk1=dvdt(v,h,n,z,t,spk);
    hk1=dhdt(v,h);
    nk1=dndt(v,n);
    zk1=dzdt(v,z);
    
    vk2=dvdt(v+vk1*dt/2,h+hk1*dt/2,n+nk1*dt/2,z+zk1*dt/2,t+dt/2,spk);
    hk2=dhdt(v+vk1*dt/2,h+hk1*dt/2);
    nk2=dndt(v+vk1*dt/2,n+nk1*dt/2);
    zk2=dzdt(v+vk1*dt/2,z+zk1*dt/2);
    
    vk3=dvdt(v+vk2*dt/2,h+hk2*dt/2,n+nk2*dt/2,z+zk2*dt/2,t+dt/2,spk);
    hk3=dhdt(v+vk2*dt/2,h+hk2*dt/2);
    nk3=dndt(v+vk2*dt/2,n+nk2*dt/2);
    zk3=dzdt(v+vk2*dt/2,z+zk2*dt/2);
    
    vk4=dvdt(v+vk3*dt,h+hk3*dt,n+nk3*dt,z+zk3*dt,t+dt,spk);
    hk4=dhdt(v+vk3*dt,h+hk3*dt);
    nk4=dndt(v+vk3*dt,n+nk3*dt);
    zk4=dzdt(v+vk3*dt,z+zk3*dt);
    
    v=v+(dt/6)*(vk1+2*vk2+2*vk3+vk4);
    h=h+(dt/6)*(hk1+2*hk2+2*hk3+hk4);
    n=n+(dt/6)*(nk1+2*nk2+2*nk3+nk4);
    z=z+(dt/6)*(zk1+2*zk2+2*zk3+zk4);
    vt3=v;

    %save the time and cell index for each spike
    for j=1:N
        if vt2(j)>=0 && vt2(j)>=vt1(j) && vt2(j)>=vt3(j)
            spktemp=t+0.001;
            spk(3,j)=spk(2,j);
            spk(2,j)=spk(1,j);
            spk(1,j)=spktemp;
            spktime(ind)=spktemp;
            spkcell(ind)=j;
            ind=ind+1;
        end
    end
end

spktime(spktime==0)=[];
spkcell(spkcell==0)=[];

%generate raster plot with sorting
[~,Index_ex]=sort(Iapp_ex,'descend');
figure(1)
for s=1:num_ex
    loc=find(spkcell==Index_ex(s));
    ti=spktime(loc);
    scatter(ti,s*ones(1,length(ti)),'.','r');
    hold on
end
for s=(num_ex+1):N
    loc=find(spkcell==s);
    ti=spktime(loc);
    scatter(ti,s*ones(1,length(ti)),'.','b');
    hold on
end

hold off

end