%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Robust Profile Clustering 
%% Programmer: BJKS     
%% Data: HCHS/SOL
%% Subpop: Center/site
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
nhanesfpedadult=readtable('nhanes_lowFadultdata1Dec2021.csv');

%remove NaN rows

SEQN=nhanesfpedadult.SEQN;
dietwt=nhanesfpedadult.dietwt8yr;
food=table2array(nhanesfpedadult(:,37:65));
age=nhanesfpedadult.RIDAGEYR;
race=nhanesfpedadult.RIDRETH3;
edu=nhanesfpedadult.DMDEDUC2;
marital=nhanesfpedadult.DMDMARTL;


n=length(SEQN);

%Create subpopulations 
subpop=NaN(n,1);
for i=1:n

    if race(i)==1
        subpop(i)=1;
    elseif race(i)==2
        subpop(i)=2;
    elseif race(i)==3
        subpop(i)=3;
    elseif race(i)==4
        subpop(i)=4;
    elseif race(i)==6
        subpop(i)=5;
    end
end

data_in=[transpose(1:n) SEQN subpop food];

data_in(any(isnan(data_in),2),:)=[];


n0=size(data_in,1);
food_hr=data_in(:,4:end);
subpop_i=data_in(:,3);
keepid=data_in(:,1);




    k_max=30;
    n=size(food_hr,1);
    p=size(food_hr,2);
    d_max=max(food_hr(:));
    d=max(food_hr);
    S=length(unique(subpop_i));
    n_s=zeros(S,1);
    for s=1:S
        n_s(s)=length(food_hr(subpop_i==s,:));
    end

    %vectorization of data
     idz = repmat(1:p,n,1); idz = idz(:);
     y_d = food_hr(:); lin_idx = sub2ind([p,d_max],idz,y_d);
        %subvectorization
        idz_s=cell(S,1);
        idz_s{S}=[];
        lin_idxS=cell(S,1);
        lin_idxS{S}=[];
        y_s=cell(S,1);
        y_s{S}=[];
     for s=1:S
         idzs=repmat(1:p,n_s(s),1); 
         idz_s{s}=idzs(:);
         food_s=food_hr(subpop_i==(s),:);
         ys=food_s(:);
         y_s{s}=ys;
         lin_idxS{s}=sub2ind([p,d_max],idz_s{s},ys);
     end

    %% SET UP PRIORS %%
   
        %beta - beta-bernoulli process 
    abe=1; bbe=1; %hypers for beta
    beta=ones(1,S);

        %nu_j^s - probability of allocation global/local
        beta_s=1;
    nu=betarnd(1,beta_s,[S,p]);

        %pi_h for all classes
a_pi=ones(1,k_max)/100;
pi_h=drchrnd(a_pi,1);

        %phi - cluster index

        rr = unifrnd(0,1,[n,1]);
       pisums=[0 cumsum(pi_h)];
       C_i=zeros(n,1);
x_ci=zeros(n,k_max);
        for l = 1:k_max
            ind = rr>pisums(l) & rr<=pisums(l+1);
            C_i(ind==1) = l;
            x_ci(:,l)=ind;
        end 
n_C_i=sum(x_ci);


          %global theta0/1
     eta=ones(1,d_max);
    theta0=zeros(p,k_max,d_max);
    theta1=zeros(S,p,k_max,d_max);

    for k=1:k_max
        for j=1:p
            dj=d(j);
            theta0(j,k,1:dj)=drchrnd(eta(1:dj),1);
            for s=1:S
            theta1(s,j,k,1:dj)=drchrnd(eta(1:dj),1);
            end
        end
    end

        %global G_ij
        G_ij=zeros(n,p);
     for s=1:S
         ns=n_s(s);
         nu_s=nu(s,:);
         G_ij(subpop_i==(s),:) = repmat(binornd(1,nu_s),[ns,1]);     % family index (0=global family,1=local family)
     end


    %% SUBPOPULATION LPP NESTS %%
    %determine number of global diets in each subpopulation
lambda_sk=drchrnd(a_pi,S);

        L_ij=zeros(n,p);
    for s=1:S
        rs = unifrnd(0,1,[n_s(s),p]); 
        ws_sums=[0 cumsum(lambda_sk(s,:))];
        phis=zeros(n_s(s),p);
        for l = 1:k_max
            ind = rs>ws_sums(l) & rs<=ws_sums(l+1);
            phis(ind==1) = l;
        end 
        L_ij(subpop_i==(s),:)=phis;
    end



    %% ------------ %%
    %% data storage %%
    %% ------------ %%
    nrun=25000;
    beta_out=zeros(nrun,S);
    nu_out=zeros(nrun,S,p);
    pi_out=zeros(nrun,k_max);
    thin=5;
    theta0_out=zeros(nrun/thin,p,k_max,d_max);
    theta1_out=zeros(nrun/thin,S,p,k_max,d_max);
    lambda_out=zeros(nrun,S,k_max);
    ci_out=zeros(nrun,n);
    n_Lij=zeros(S,k_max);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% POSTERIOR COMPUTATION %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    As=zeros(n,p);
    p_ij=zeros(n,p);
tic
    for iter=1:nrun

        %% -- update z prob -- %%
     for s=1:S
         phrow1=L_ij(subpop_i==(s),:); phrow1=phrow1(:);
         theta1s=reshape(theta1(s,:,:,:),[p,k_max,d_max]);
         A = theta1s(sub2ind([p,k_max,d_max],idz_s{s},phrow1,y_s{s})); %index of subjects in theta1class
        A = reshape(A,[n_s(s),p]);
        As(subpop_i==(s),:)=A;
     end

     phrow0=repmat(C_i,[1,p]); phrow0=phrow0(:);
     B = theta0(sub2ind([p,k_max,d_max],idz,phrow0,y_d));
     B = reshape(B,[n,p]);

     for s=1:S
         ns=n_s(s); nu_s=nu(s,:);
         p_ij(subpop_i==(s),:)=(repmat(nu_s,[ns,1]).*B(subpop_i==(s),:)./((repmat(nu_s,[ns,1]).*B(subpop_i==(s),:))+(repmat(1-nu_s,[ns,1]).*As(subpop_i==(s),:))));
     end

        G_ij=binornd(1,p_ij);

        %% -- update pi_h -- %%
         for h=1:k_max
             n_C_i(h)=sum(C_i==h);
         end

        a_pih=a_pi+n_C_i;
        pi_h=drchrnd(a_pih,1);

        pi_out(iter,1:length(pi_h))=pi_h;

      %%-- update lambda_sk --%%

  for s=1:S
      L_s=L_ij(subpop_i==(s),:);

      for l=1:k_max
        n_Lij(s,l)=sum(L_s(:)==l);
%          n_phij(s,l)=sum(phi_sz(:)==l);
      end
      kn_Lij=a_pi+n_Lij(s,:);
      lambda_sk(s,:)=drchrnd(kn_Lij,1);
  end

        lambda_out(iter,:,:)=lambda_sk;



      %% -- phi ~multinomial(pi_h) -- %%

  Cp_k=zeros(n,k_max);
for k=1:k_max
    t0h=reshape(theta0(:,k,:),p,d_max);
    tmpmat0=reshape(t0h(lin_idx),[n,p]);
    Cp_k(:,k)=pi_h(k)*prod(tmpmat0.^G_ij,2);
end
probCi = bsxfun(@times,Cp_k,1./(sum(Cp_k,2)));
    x_ci=mnrnd(1,probCi); [r, c]=find(x_ci); x_gc=[r c];
    x_gc=sortrows(x_gc,1); C_i=x_gc(:,2);
for s=1:S
    Lijs=zeros(n_s(s),k_max,p);
         for h = 1:k_max
            theta1hs = reshape(theta1(s,:,h,:),p,d_max);
            tmpmat1 = reshape(theta1hs(lin_idxS{s}),[n_s(s),p]);
             Lijs(:,h,:) = lambda_sk(s,h) * tmpmat1.^(G_ij(subpop_i==(s),:)==0);
         end  
        sumLijs=repmat(sum(Lijs,2),[1,k_max,1]);
        zupS = Lijs./sumLijs;
        for j=1:p
            sub_pj=reshape(zupS(:,:,j),[n_s(s),k_max]);
            L_sub=mnrnd(1,sub_pj);
            [r, c]=find(L_sub); x_l=[r c];
            x_sl=sortrows(x_l,1);
            L_ij(subpop_i==(s),j)=x_sl(:,2);
        end
end
%store global cluster index for postprocess/relabelling
ci_out(iter,:)=C_i;

        % - update theta - %
    dmat0=zeros(p,d_max);
    dmat1=zeros(p,d_max);
    for k=1:k_max
        C_is=repmat(C_i,[1,p]).*G_ij;
         ph0 = (C_is==k); %subj's in global cluster h
            for c = 1:d_max
                 dmat0(:,c) = sum((food_hr==c).*ph0)';
            end
            for j=1:p
                dj=d(j);
                a_tn0=eta(1:dj)+dmat0(j,1:dj);
                theta0(j,k,1:dj) = drchrnd(a_tn0,1);
            end
    end

      for s=1:S
          phis=L_ij(subpop_i==(s),:).*(1-G_ij(subpop_i==(s),:));
          foods=food_hr(subpop_i==(s),:);
          for l=1:k_max
              ph1=(phis==l);
              for c=1:d_max
                dmat1(:,c) = sum((foods==c).*ph1);
              end
            for j=1:p
                dj=d(j);
                a_tn1=eta(1:dj)+dmat1(j,1:dj);
                theta1(s,j,l,1:dj) = drchrnd(a_tn1,1);
            end 
          end
      end
      if mod(iter,thin)==0
         theta0_out(iter/thin,:,1:size(theta0,2),:)=theta0;
         theta1_out(iter/thin,:,:,1:size(theta1,3),:)=theta1;
      end



        % update nu_j %
        for s=1:S
            Gs=G_ij(subpop_i==(s),:);
            nu(s,:) = betarnd(1 + sum(Gs), beta(s) + sum(1-Gs));
        end
      nu(nu==1) = 1-1e-06;
      nu(nu==0) = 1e-06; 

     nu_out(iter,:,:)=nu;

      % - update beta - %
      for s=1:S,
        beta(s) = gamrnd(abe + p,1./( bbe - sum(log(1-nu(s,:)))));
      end
     beta_out(iter,:)=beta;

%% RELABELLING STEP TO ENCOURAGE MIXING %%
if mod(iter,10)==0
new_order=randperm(k_max);
newC_i=C_i;
for k=1:k_max
newC_i(C_i==k)=new_order(k);
end
C_i=newC_i;
theta0=theta0(:,new_order,:);
end




    end
eltime=toc;
    burn=nrun/5;
    beta_burn=beta_out(burn+1:end,:);
    lambda_burn=lambda_out(burn+1:end,:,:);
    pi_burn=pi_out(burn+1:end,:);
    theta0_burn=theta0_out((burn/thin)+1:end,:,:,:);
    theta1_burn=theta1_out((burn/thin)+1:end,:,:,:,:);
    nu_burn=nu_out(burn+1:end,:,:);
    ci_burn=ci_out(burn+1:end,:);

    %posterior median of beta
    beta_med=median(beta_burn);
   
    nu=reshape(median(nu_burn,1),[S,p]);

%posterior medians of local cluster model
lambdas=reshape(median(lambda_burn,1),[S,k_max]);
theta1=reshape(median(theta1_burn,1),[S,p,k_max,d_max]);

lambdas_x=cell(S,1); lambdas_x{S}=[];
theta1_x=cell(S,1); theta1_x{S}=[];
for s=1:S
    wss=lambdas(s,:);
lambdas_x{s}=wss(wss>0.05);
theta1_x{s}=reshape(theta1(s,:,wss>0.05,:),[p,length(lambdas_x{s}),d_max]);

end

nu_med=reshape(median(nu_burn),[S,p]);

save('RPCnhanesLOWF_MCMCout','beta_burn','lambdas_x','pi_burn','ci_burn','theta0_burn','theta1_burn','nu_burn','eltime','-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESS: PASAPILIOPOULIS SWITCH %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=size(pi_burn,1); m_perm=size(theta0_burn,1); m_thin=m/m_perm;
k_med=median(sum(pi_burn>0.05,2));
pd=pdist(transpose(ci_burn),'hamming'); %prcnt differ
cdiff=squareform(pd); %Dij
Zci=linkage(cdiff,'complete');
figure; %save dendrogram of hierarchical clustering
dendrogram(Zci);
saveas(gcf,'nhanesRPC_dendrogram.png')

clust0 = cluster(Zci,'maxclust',k_med); %choose k0=5;

%Ordered MCMC set
ci_relabel=zeros(m,k_med); %ci_relabelin=zeros(m,k_in); 
for l=1:k_med
    ci_relabel(:,l) = mode(ci_burn(:,clust0==l),2);
end

pi_order=zeros(m,k_med);
theta0_order=zeros(m_perm,p,k_med,d_max);
for iter=1:m
    iter_order=ci_relabel(iter,:);
    pi_h1=pi_burn(iter,iter_order);
    pi_order(iter,:)=pi_h1/sum(pi_h1);
    if mod(iter,m_thin)==0
        iter_thin=iter/m_thin;
        theta0_order(iter_thin,:,:,:)=theta0_burn(iter_thin,:,iter_order,:);
    end
end

theta0_med=reshape(median(theta0_order),[p,k_med,d_max]);
pi_med=median(pi_order);
pi_med=pi_med/sum(pi_med);

save('RPCnhanesLOWF_MCMCmed','beta_med','lambdas_x','pi_med','theta0_med','theta1_x','nu_med','eltime','Zci','-v7.3');

[p,k_med,d]=size(theta0_med);
ks=zeros(S,1);
val1=cell(S,1); ind1=cell(S,1);
theta1med=cell(S,1);
theta1_mat=zeros(p,S,d_max);
for s=1:S
    ks(s)=size(theta1_x{s},2);
    theta1s=theta1_x{s};
    [val1{s},ind1{s}]=max(theta1s,[],3);
    theta1s=theta1s./sum(theta1s,3);
    theta1med{s}=reshape(theta1s,[p,d_max]);
    %theta1_x{s}=theta1_x{s}./sum(theta1_x{s});
    theta1_mat(:,s,:)=theta1med{s};
end


%%Modal pattern of Global Clusters %%
theta0_medi=theta0_med./sum(theta0_med,3);

[val0,ind0]=max(theta0_medi,[],3);

theta0_probs=cell(d_max,1);
for dk=1:d_max
    theta0_probs{dk}=reshape(theta0_med(:,:,dk),[p,k_med]);
end


%%

flabels={'Citrus, Melon Berries','Other fruit','Fruit juice','Dk Green Veg',...
    'Tomatoes','Other Red/Org veg','Potatoes','Other starchy veg','Other veg',...
    'Legumes (veg)','Whole grains','Refined grains','Meat(ns)','Cured meats',...
    'Organ meat','Poultry','Seafood (highn3)','Seafood (lown3)','Eggs','Soybean',...
    'Nuts/seeds','Legumes (protein)','Milk','Yogurt','Cheese','Oils','Solid fat',...
    'Added sugar','Alcohol'};

%stacked and grouped bar plot for foods 
figure; 
plotBarStackGroupsLab(theta1_mat,flabels)
cb = colorbar;
set(cb,'Limits', [0,1], 'Ticks',.125:.25:1,'TickLabels',["None" "Low" "Medium" "High"])

%split figure for foods;
theta1mat_1=theta1_mat(1:15,:,:);
flabels_1=flabels(1:15); 
plotBarStackGroupsLab(theta1mat_1,flabels_1)
cb = colorbar;
set(cb,'Limits', [0,1], 'Ticks',.125:.25:1,'TickLabels',["None" "Low" "Medium" "High"])

theta1mat_2=theta1_mat(16:29,:,:);
flabels_2=flabels(16:29);
plotBarStackGroupsLab(theta1mat_2,flabels_2)
cb = colorbar;
set(cb,'Limits', [0,1], 'Ticks',.125:.25:1,'TickLabels',["None" "Low" "Medium" "High"])



% Posterior probability plot of no consumption
figure; bar(transpose(theta1_mat(:,:,1)))
xticks(1:29)
xticklabels(flabels)
ylabel('Posterior Probability of No consumption')
saveas(gcf,'noconsum_lowadultpat.fig')


% Posterior probability plot of high consumption
figure; plot(theta0_probs{4},'Linewidth',1)
xticks(1:29)
xticklabels(flabels)
yticks(0:0.1:1)
ylabel('Posterior Probability of High consumption')
saveas(gcf,'highconsum_lowadultpat.fig')

figure; 
    h=heatmap(ind0)
    h.YDisplayLabels = flabels;
    h.XLabel = "Dietary Profile";
    h.Colormap = parula
saveas(gcf,'theta0_lowFpattern.fig')
    
format long
nu_t=transpose(round(nu_med*100,1));
%Plot of global/local allocation probability 
racelab = {'Mexican', 'Other Hispanic', 'NH White', 'NH Black', 'NH Asian'}


for s=1:S
    figure; barh(theta1med{s},'stacked')
    xlim(0:1)
    yticks(1:p)
    yticklabels(flabels)
    title(strcat('Dietary Profile - ',racelab{s}))
    xlabel('Posterior Probability')
    legend({'None','Low','Med','High'},'Location','southwestoutside')
end



figure; 
    hn=heatmap(nu_t)
    hn.YDisplayLabels = flabels;
    hn.XDisplayLabels = racelab;
saveas(gcf,'nu_lowwomen.fig')


%% save final assignments

G_med=zeros(n,p);
    for s=1:S
       G_med(subpop_i==s,:)=binornd(1,repmat(nu_med(s,:),[n_s(s),1]));
    end
    
delmed=zeros(n,k_med);
   for h = 1:k_med
        t0h = reshape(theta0_med(:,h,:),[p,d_max]);
        theta0h=bsxfun(@times,t0h,1./sum(t0h,2));
        tmpmat0 = reshape(theta0h(lin_idx),[n,p]);
        delmed(:,h) = pi_med(h)*prod(tmpmat0.^G_med,2);
    end 
    zup0 = bsxfun(@times,delmed,1./(sum(delmed,2)));
    [p_ci,rpc_ci]=max(zup0,[],2);

seqn_rpc=SEQN(keepid);
dietwt8yr=dietwt(keepid);
nhanes_rpc=table(seqn_rpc,dietwt8yr,subpop_i,rpc_ci);
writetable(nhanes_rpc,'NHANESLowFrpc_assign1Dec2021.csv');
