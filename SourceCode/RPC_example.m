%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Robust Profile Clustering Method 
%Programmer: Briana Stephenson
%Data: Simulated example 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;
    %unload data
    load('ExampleData');
    dataf=subdata;
    subpop=subpop_i;

    k_max=50;
    n=size(dataf,1);
    p=size(dataf,2);
    d=max(dataf(:));
    S=length(unique(subpop));
    n_s=zeros(S,1);
    for s=1:S
        n_s(s)=length(dataf(subpop==s,:));
    end

    %vectorization of data
     idz = repmat(1:p,n,1); idz = idz(:);
     y_d = dataf(:); lin_idx = sub2ind([p,d],idz,y_d);
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
         food_s=dataf(subpop==(s),:);
         ys=food_s(:);
         y_s{s}=ys;
         lin_idxS{s}=sub2ind([p,d],idz_s{s},ys);
     end

    %% SET UP PRIORS %%
   
        %beta - beta-bernoulli process 
    a_beta=1; b_beta=1; %hypers for beta
    beta=ones(1,S);

        %nu_j^s - probability of allocation global/local
    beta_s=1;
    nu=betarnd(1,beta_s,[S,p]);

        %pi_h for all classes
a_pi=ones(1,k_max)/k_max;
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
     eta=ones(1,d);
    theta0=zeros(p,k_max,d);
    theta1=zeros(S,p,k_max,d);

    for k=1:k_max,
        for j=1:p,
            theta0(j,k,:)=drchrnd(eta,1);
            for s=1:S,
            theta1(s,j,k,:)=drchrnd(eta,1);
            end
        end
    end

        %global G_ij
        G_ij=zeros(n,p);
     for s=1:S
         ns=n_s(s);
         nu_s=nu(s,:);
         G_ij(subpop==(s),:) = repmat(binornd(1,nu_s),[ns,1]);     % family index (0=global family,1=local family)
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
        L_ij(subpop==(s),:)=phis;
    end



    %% ------------ %%
    %% data storage %%
    %% ------------ %%
    nrun=25000;
    beta_out=zeros(nrun,S);
    nu_out=zeros(nrun,S,p);
    pi_out=zeros(nrun,k_max);
    thin=100;
    theta0_out=zeros(nrun/thin,p,k_max,d);
    theta1_out=zeros(nrun/thin,S,p,k_max,d);
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
         phrow1=L_ij(subpop==(s),:); phrow1=phrow1(:);
         theta1s=reshape(theta1(s,:,:,:),[p,k_max,d]);
         A = theta1s(sub2ind([p,k_max,d],idz_s{s},phrow1,y_s{s})); %index of subjects in theta1class
        A = reshape(A,[n_s(s),p]);
        As(subpop==(s),:)=A;
     end

     phrow0=repmat(C_i,[1,p]); phrow0=phrow0(:);
     B = theta0(sub2ind([p,k_max,d],idz,phrow0,y_d));
     B = reshape(B,[n,p]);

     for s=1:S
         ns=n_s(s); nu_s=nu(s,:);
         p_ij(subpop==(s),:)=(repmat(nu_s,[ns,1]).*B(subpop==(s),:)./((repmat(nu_s,[ns,1]).*B(subpop==(s),:))+(repmat(1-nu_s,[ns,1]).*As(subpop==(s),:))));
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
      L_s=L_ij(subpop==(s),:);

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
    t0h=reshape(theta0(:,k,:),p,d);
    tmpmat0=reshape(t0h(lin_idx),[n,p]);
    Cp_k(:,k)=pi_h(k)*prod(tmpmat0.^G_ij,2);
end
probCi = bsxfun(@times,Cp_k,1./(sum(Cp_k,2)));
    x_ci=mnrnd(1,probCi); [r, c]=find(x_ci); x_gc=[r c];
    x_gc=sortrows(x_gc,1); C_i=x_gc(:,2);
for s=1:S
    Lijs=zeros(n_s(s),k_max,p);
         for h = 1:k_max
            theta1hs = reshape(theta1(s,:,h,:),p,d);
            tmpmat1 = reshape(theta1hs(lin_idxS{s}),[n_s(s),p]);
             Lijs(:,h,:) = lambda_sk(s,h) * tmpmat1.^(G_ij(subpop==(s),:)==0);
         end  
        sumLijs=repmat(sum(Lijs,2),[1,k_max,1]);
        zupS = Lijs./sumLijs;
        for j=1:p
            sub_pj=reshape(zupS(:,:,j),[n_s(s),k_max]);
            L_sub=mnrnd(1,sub_pj);
            [r, c]=find(L_sub); x_l=[r c];
            x_sl=sortrows(x_l,1);
            L_ij(subpop==(s),j)=x_sl(:,2);
        end
end
%store global cluster index for postprocess/relabelling
ci_out(iter,:)=C_i;

        % - update theta - %
    dmat0=zeros(p,d);
    dmat1=zeros(p,d);
    for k=1:k_max
        C_is=repmat(C_i,[1,p]).*G_ij;
         ph0 = (C_is==k); %subj's in global cluster h
            for c = 1:d
                 dmat0(:,c) = sum((dataf==c).*ph0)';
            end
            for j=1:p
                a_tn0=eta+dmat0(j,:);
                theta0(j,k,:) = drchrnd(a_tn0,1);
            end
    end

      for s=1:S
          phis=L_ij(subpop==(s),:).*(1-G_ij(subpop==(s),:));
          foods=dataf(subpop==(s),:);
          for l=1:L_s(s)
              ph1=(phis==l);
              for c=1:d
                dmat1(:,c) = sum((foods==c).*ph1);
              end
            for j=1:p
                a_tn1=eta+dmat1(j,:);
                theta1(s,j,l,:) = drchrnd(a_tn1,1);
            end 
          end
      end
      if mod(iter,thin)==0
         theta0_out(iter/thin,:,1:size(theta0,2),:)=theta0;
         theta1_out(iter/thin,:,:,1:size(theta1,3),:)=theta1;
      end



        % update nu_j %
        for s=1:S
            Gs=G_ij(subpop==(s),:);
            nu(s,:) = betarnd(1 + sum(Gs), beta(s) + sum(1-Gs));
        end
      nu(nu==1) = 1-1e-06;
      nu(nu==0) = 1e-06; 

     nu_out(iter,:,:)=nu;

      % - update beta - %
      for s=1:S
        beta(s) = gamrnd(a_beta + p,1./( b_beta - sum(log(1-nu(s,:)))));
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
   
        figure; %check mixing of pi parameter
        plot(pi_burn)
        saveas(gcf,'pis.png')
        figure; %check mixing of beta parameter
        plot(beta_burn)
        saveas(gcf,'betas.png')

    beta=median(beta_burn);
    nu=reshape(median(nu_burn,1),[S,p]);
    
    clf
subplot(1,2,1);
heatmap(transpose(nu), 1:S,1:p);
title('Derived Local Deviations');

subplot(1,2,2);
heatmap(sub_nu, 1:S, 1:p);
title('True Local Deviations');
saveas(gcf,strcat('G_deviations.png'))

%posterior medians of local cluster model
lambdas=reshape(median(lambda_burn,1),[S,k_max]);
theta1=reshape(median(theta1_burn,1),[S,p,k_max,d]);

lambdas_x=cell(S,1); lambdas_x{S}=[];
theta1_x=cell(S,1); theta1_x{S}=[];
for s=1:S
    wss=lambdas(s,:);
lambdas_x{s}=wss(wss>0.01);
theta1_x{s}=reshape(theta1(s,:,wss>0.01,:),[p,length(lambdas_x{s}),d]);
end

nu_med=reshape(median(nu_burn),[S,p]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESS: PASAPILIOPOULIS SWITCH %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=size(ci_burn,1);
k_med=median(sum(pi_burn>0.05,2));
pd=pdist(transpose(ci_burn),'hamming'); %prcnt differ
cdiff=squareform(pd); %Dij
Zci=linkage(cdiff,'complete');
figure; %save dendrogram of hierarchical clustering
dendrogram(Zci);
saveas(gcf,'dendrogram.png')

clust_med = cluster(Zci,'maxclust',k_med);

%Ordered MCMC set
ci_relabel=zeros(m,k_med); %ci_relabelin=zeros(m,k_in); 
for l=1:k_med
    ci_relabel(:,l) = mode(ci_burn(:,clust_med==l),2);
end
m_perm=m/thin;
pi_order=zeros(m,k_med);
theta0_order=zeros(m_perm,p,k_med,d);
for iter=1:m
    iter_order=ci_relabel(iter,:);
    pi_h1=pi_burn(iter,iter_order);
    pi_order(iter,:)=pi_h1/sum(pi_h1);
    if mod(iter,thin)==0
        iter_thin=iter/thin;
        theta0_order(iter_thin,:,:,:)=theta0_burn(iter_thin,:,iter_order,:);
    end
end

theta0_med=reshape(median(theta0_order),[p,k_med,d]);
pi_med=median(pi_order);
pi_med=pi_med/sum(pi_med);

save('RPCparm_ex','beta_med','lambdas_x','pi_med','theta0_med','theta1_x','nu_med','eltime','-v7.3');

%%EXPORT VARIABLS TO EXCEL FILE
%filename = 'RPCparameters.xlsx';
dlmwrite('betas.txt',beta_med);
dlmwrite('pis.txt',pi_med);
nu_t=transpose(nu_med);
dlmwrite('nus.txt',nu_t);

%create modal pattern for global clusters
[global_val, global_patt]=max(theta0_med,[],3);
dlmwrite('GlobalPattern_predicted.txt',global_patt);
