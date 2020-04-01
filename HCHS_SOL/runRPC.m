%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Robust Profile Clustering 
%% Programmer: BJKS     
%% Data: HCHS/SOL
%% Subpop: Center/site
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
fpq_data=xlsread('HCHSrpc_exampledata.xlsx');

food=fpq_data(:,3:end);
subpop=fpq_data(:,2);
id=fpq_data(:,1);


    k_max=50;
    n=size(food,1);
    p=size(food,2);
    d_max=max(food(:));
    d=max(food);
    S=length(unique(subpop));
    n_s=zeros(S,1);
    for s=1:S
        n_s(s)=sum(subpop==s);
    end

    %vectorization of data
     idz = repmat(1:p,n,1); idz = idz(:);
     y_d = food(:); lin_idx = sub2ind([p,d_max],idz,y_d);
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
         food_s=food(subpop==(s),:);
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
     eta=ones(1,d_max);
    theta0=zeros(p,k_max,d_max);
    theta1_med=zeros(S,p,k_max,d_max);

    for k=1:k_max
        for j=1:p
            dj=d(j);
            theta0(j,k,1:dj)=drchrnd(eta(1:dj),1);
            for s=1:S
            theta1_med(s,j,k,1:dj)=drchrnd(eta(1:dj),1);
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


    %% SUBPOPULATION PROFILE %%
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
    nrun=2500;
    beta_out=zeros(nrun,S);
    nu_out=zeros(nrun,S,p);
    pi_out=zeros(nrun,k_max);
    thin=25;
    theta0_out=zeros(nrun/thin,p,k_max,d_max);
    theta1_out=zeros(nrun/thin,S,p,k_max,d_max);
    lambda_out=zeros(nrun,S,k_max);
    ci_out=zeros(nrun,n);
    n_Lij=zeros(S,k_max);
    log_locrpc=zeros(n,1);
    dic_out=zeros(n,1);

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
         theta1s=reshape(theta1_med(s,:,:,:),[p,k_max,d_max]);
         A = theta1s(sub2ind([p,k_max,d_max],idz_s{s},phrow1,y_s{s})); %index of subjects in theta1class
        A = reshape(A,[n_s(s),p]);
        As(subpop==(s),:)=A;
     end

     phrow0=repmat(C_i,[1,p]); phrow0=phrow0(:);
     B = theta0(sub2ind([p,k_max,d_max],idz,phrow0,y_d));
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
    t0h=reshape(theta0(:,k,:),p,d_max);
    tmpmat0=reshape(t0h(lin_idx),[n,p]);
    Cp_k(:,k)=pi_h(k)*prod(tmpmat0.^G_ij,2);
end
    log_globrpc=sum(log(Cp_k),2);
    probCi = bsxfun(@times,Cp_k,1./(sum(Cp_k,2)));
    x_ci=mnrnd(1,probCi); [r, c]=find(x_ci); x_gc=[r c];
    x_gc=sortrows(x_gc,1); C_i=x_gc(:,2);
for s=1:S
    Lijs=zeros(n_s(s),k_max,p);
         for h = 1:k_max
            theta1hs = reshape(theta1_med(s,:,h,:),p,d_max);
            tmpmat1 = reshape(theta1hs(lin_idxS{s}),[n_s(s),p]);
             Lijs(:,h,:) = lambda_sk(s,h) * tmpmat1.^(G_ij(subpop==(s),:)==0);
         end  
         log_locrpc(subpop==s,1)=sum(log(sum(Lijs,2)),3);
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
dic_out(iter)=sum(log_globrpc+log_locrpc);

        % - update theta - %
    dmat0=zeros(p,d_max);
    dmat1=zeros(p,d_max);
    for k=1:k_max
        C_is=repmat(C_i,[1,p]).*G_ij;
         ph0 = (C_is==k); %subj's in global cluster h
            for c = 1:d_max
                 dmat0(:,c) = sum((food==c).*ph0)';
            end
            for j=1:p
                dj=d(j);
                a_tn0=eta(1:dj)+dmat0(j,1:dj);
                theta0(j,k,1:dj) = drchrnd(a_tn0,1);
            end
    end

      for s=1:S
          phis=L_ij(subpop==(s),:).*(1-G_ij(subpop==(s),:));
          foods=food(subpop==(s),:);
          for l=1:k_max
              ph1=(phis==l);
              for c=1:d_max
                dmat1(:,c) = sum((foods==c).*ph1);
              end
            for j=1:p
                dj=d(j);
                a_tn1=eta(1:dj)+dmat1(j,1:dj);
                theta1_med(s,j,l,1:dj) = drchrnd(a_tn1,1);
            end 
          end
      end
      if mod(iter,thin)==0
         theta0_out(iter/thin,:,1:size(theta0,2),:)=theta0;
         theta1_out(iter/thin,:,:,1:size(theta1_med,3),:)=theta1_med;
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
eltime=toc/(60);
    burn=nrun/5;
    beta_burn=beta_out(burn+1:end,:);
    lambda_burn=lambda_out(burn+1:end,:,:);
    pi_burn=pi_out(burn+1:end,:);
    theta0_burn=theta0_out((burn/thin)+1:end,:,:,:);
    theta1_burn=theta1_out((burn/thin)+1:end,:,:,:,:);
    nu_burn=nu_out(burn+1:end,:,:);
    ci_burn=ci_out(burn+1:end,:);
    dic_burn=dic_out(burn+1:end);

    %posterior median of beta
    beta_med=median(beta_burn);   

%posterior medians of local cluster model
lambda_med=reshape(median(lambda_burn,1),[S,k_max]);
theta1_med=reshape(median(theta1_burn,1),[S,p,k_max,d_max]);

lambdas_x=cell(S,1); lambdas_x{S}=[];
theta1_x=cell(S,1); theta1_x{S}=[];
for s=1:S
    wss=lambda_med(s,:);
lambdas_x{s}=wss(wss>0.05);
theta1_x{s}=reshape(theta1_med(s,:,wss>0.05,:),[p,length(lambdas_x{s}),d_max]);
end

nu_t=reshape(median(nu_burn),[S,p]);

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
saveas(gcf,'hchsex_dendrogram.png')

clust0 = cluster(Zci,'maxclust',k_med); %choose k0=5;

%Reorder parameters from MCMC set
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

    % Recalculate DIC with ordered/reduced parameters
loglike_rpclocal=zeros(n,1);
dic_thin=zeros(m_perm,1);
for ii=1:m_perm
    i_thin=ii*m_thin;
    pi_thin=pi_order(i_thin,:);
    theta0_thin=reshape(theta0_order(ii,:,:,:),[p,k_med,d_max]);
    theta1_thin=reshape(theta1_burn(ii,:,:,:,:),[S,p,k_max,d_max]);
    nu_thin=reshape(nu_burn(i_thin,:,:),[S,p]);
    lambda_thin=reshape(lambda_burn(i_thin,:,:),[S,k_max]);
    pi_h=pi_thin/sum(pi_thin);
    for s=1:S
        nu_s=reshape(nu_thin(s,:),[1,p]);
        G_ij(subpop==s,:)=binornd(1,repmat(nu_s,[n_s(s),1]));
    end
    

        %assign global cluster
      Cp_k=zeros(n,k_med);
    for k=1:k_med
        t0h=reshape(theta0_thin(:,k,:),p,d_max);
        tmpmat0=reshape(t0h(lin_idx),[n,p]);
        Cp_k(:,k)=pi_h(k)*prod(tmpmat0.^G_ij,2);
    end
    probCi = bsxfun(@times,Cp_k,1./(sum(Cp_k,2)));
    log_globrpc=log(sum(Cp_k,2));
    tij_logglobrpci=sum(probCi.*log(Cp_k),2);
    w_ci=mnrnd(1,probCi); 
   for s=1:S
        lts=lambda_thin(s,:);
        lambda_ii=lts(lts>0.05);
        lambda_is=lambda_ii/sum(lambda_ii);
        ks=length(lambda_is);
        theta1_is=reshape(theta1_thin(s,:,lts>0.05,:),[p,ks,d_max]);
    Lijs=zeros(n_s(s),ks,p);
         for h = 1:ks
            theta1hs = reshape(theta1_is(:,h,:),p,d_max);
            tmpmat1 = reshape(theta1hs(lin_idxS{s}),[n_s(s),p]);
             Lijs(:,h,:) = lambda_is(h) * tmpmat1.^(G_ij(subpop==s,:)==0);
         end  
        sumLijs=repmat(sum(Lijs,2),[1,ks,1]);
        zupS = Lijs./sumLijs;
        for j=1:p
            sub_pj=reshape(zupS(:,:,j),[n_s(s),ks]);
            l_ij=mnrnd(1,sub_pj);
            [r, c]=find(l_ij); x_l=[r c];
            x_sl=sortrows(x_l,1);
            L_ij(subpop==s,j)=x_sl(:,2);
        end  
   end
    locallike=sum(Lijs,2);
 
    loglike_rpc=tij_logglobrpci; %+loglike_rpclocal;
    dic_thin(ii)=sum(loglike_rpc);

end


theta0_med=reshape(median(theta0_order),[p,k_med,d_max]);
pi_med=median(pi_order);
nu_med=reshape(median(nu_burn),[S,p]);
G_med=zeros(n,p);
    for s=1:S
       G_med(subpop==s,:)=binornd(1,repmat(nu_med(s,:),[n_s(s),1]));
    end
 like_srpcmed=zeros(n,1);
delmed=zeros(n,k_med);
   for h = 1:k_med
        t0h = reshape(theta0_med(:,h,:),[p,d_max]);
        theta0h=bsxfun(@times,t0h,1./sum(t0h,2));
        tmpmat0 = reshape(theta0h(lin_idx),[n,p]);
        delmed(:,h) = pi_med(h)*prod(tmpmat0.^G_med,2);
    end 
    zup0 = bsxfun(@times,delmed,1./(sum(delmed,2)));
    tij_logglobrpcmed=sum(zup0.*log(delmed),2);

    wm_ci=mnrnd(1,zup0); [r, c]=find(wm_ci); x_gc=[r c];
    x_gc=sortrows(x_gc,1); Ci=x_gc(:,2);
loglike_locsrpcmed=zeros(n,1);
for s=1:S
        lambda_m=lambdas_x{s};
        lambda_med=lambda_m/sum(lambda_m);
        theta1_ms=theta1_x{s};
        ks=length(lambda_med);
        
    Lijs=zeros(n_s(s),ks,p);
         for h = 1:ks
            theta1hs = reshape(theta1_ms(:,h,:),[p,d_max]);
            tmpmat1 = reshape(theta1hs(lin_idxS{s}),[n_s(s),p]);
             Lijs(:,h,:) = lambda_med(h) * tmpmat1.^(G_med(subpop==s,:)==0);
         end  
        sumLijs=repmat(sum(Lijs,2),[1,ks,1]);
        zupS = Lijs./sumLijs;
        for j=1:p
            sub_pj=reshape(zupS(:,:,j),[n_s(s),ks]);
            l_ij=mnrnd(1,sub_pj);
            [r, c]=find(l_ij); x_l=[r c];
            x_sl=sortrows(x_l,1);
            L_ij(subpop==s,j)=x_sl(:,2);
        end  
    locallike_med=sum(Lijs,2);
    loglike_locsrpcmed(subpop==s)=sum(log(locallike_med),3);
end

dic_med=tij_logglobrpcmed+loglike_locsrpcmed;


dic=-4*median(dic_thin)+2*sum(dic_med);


%% modal patterns of local profile %%
[p,k_med,d]=size(theta0_med);
ks=zeros(S,1);
theta1val=cell(S,1); theta1ind=cell(S,1);
for s=1:S
    ks(s)=size(theta1_x{s},2);
    [theta1val{s},theta1ind{s}]=max(theta1_x{s},[],3);
end


%%Modal pattern of Global Clusters %%
t0med=bsxfun(@times,theta0_med,1./sum(theta0_med,3));
[theta0val,theta0ind]=max(t0med,[],3);


save('RPCexample_MCMCoutput','beta_med','lambdas_x','pi_med','theta0val','theta0ind','theta1val','theta1ind','theta1_x','nu_med','Zci','-v7.3');

%assign global profile based on posterior medians
[rpc_p,rpc_i]=max(zup0,[],2);

%%generate output table%%

hchs_rpc=table(id,subpop,rpc_i,rpc_p);
writetable(hchs_rpc,'HCHSexample_RPCoutput.xlsx');
