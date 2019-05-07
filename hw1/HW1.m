
function HW1()
    clear;
    clc;
    
    vol=0.1*sqrt(1/260);
    lmbd=60*60*24;
    dlt_t=0.1/lmbd;
    trdng_prob=1-exp(-lmbd*dlt_t);
    sprd_clnt=1e-4;
    sprd_dlr=2e-4;
    dlt_lmt=3;
    
    n_stps=500;
    n_rns=1e4;
    
    run_sim(vol,sprd_clnt,sprd_dlr,dlt_lmt,dlt_t,n_stps,n_rns,trdng_prob,0);
    run_sim(vol,sprd_clnt,sprd_dlr,dlt_lmt,dlt_t,n_stps,n_rns,trdng_prob,1);

end

function res=run_sim(vol,sprd_clnt,sprd_dlr,dlt_lmt,dlt_t,n_stps,n_rns,trdng_prob,hdg_type)

    tic
    rng(1);
    pnls=NaN*ones(1,n_rns);
    trades=NaN*ones(1,n_rns);
    hedges=NaN*ones(1,n_rns);
    
    for this_run=1:n_rns

       hdg_bool=zeros(n_stps,1);
       trd_bool=zeros(n_stps,1);
       
       nrml_rnd=random('Normal',0,sqrt(dlt_t),[n_stps,1]);
       trd_rnd=random('Uniform',0,1,[n_stps,1]);
       pos_rnd=random('Binomial',1,0.5,[n_stps,1]);
       pos_rnd(pos_rnd==0)=-1;
       
       trd_bool(trd_rnd<trdng_prob)=1;
       pos_rnd(~trd_bool)=0;
       postn=cumsum(pos_rnd);
       hdg_wht=[];
       
       while ~isempty(cond_hedge(postn,dlt_lmt,hdg_type))
           hdg_this=cond_hedge(postn,dlt_lmt,hdg_type);
           hdg_bool(hdg_this)=1;
           sgn=postn(hdg_this)/abs(postn(hdg_this));
           hdg_wht=[hdg_wht postn(hdg_this)];
           postn(hdg_this:end)=cumsum([sgn*dlt_lmt*hdg_type;pos_rnd(hdg_this+1:end)]);
       end
       sgn=(postn+eps)./abs(postn+eps);
       postn(hdg_bool==1)=hdg_wht;
     
       hdg_ntnl=hdg_type*sgn.*postn+((-1)^hdg_type)*dlt_lmt;
       
       pnl_hdges=-hdg_bool.*hdg_ntnl.*sprd_dlr.*0.5.*cumprod(1+vol*nrml_rnd(:));
       pnl_trdes=trd_bool.*sprd_clnt*0.5.*cumprod(1+vol*nrml_rnd(:));
       pnl_iters=postn.*vol.*cumprod(1+vol*nrml_rnd(:)).*nrml_rnd;
       
       pnl=cumsum(pnl_hdges+pnl_trdes+pnl_iters);
       pnls(this_run)=pnl(end);
       trades(this_run)=sum(trd_bool);
       hedges(this_run)=sum(hdg_bool);
        
    end
    
    res.nruns=n_rns;
    res.sharpe=mean(pnls)/std(pnls);
    res.mean_pnls=mean(pnls);
    res.std_pnls=std(pnls);
    res.mean_trades=mean(trades);
    res.mean_hedges=mean(hedges);
    
    toc;res
end

function hedge_this=cond_hedge(postn,dlt_lmt,hdg_type)
    if hdg_type==0
        hedge_this=find(abs(postn)>=dlt_lmt,1);
    elseif hdg_type==1
        hedge_this=find(abs(postn)>dlt_lmt,1);
    else
        disp();
    end
end