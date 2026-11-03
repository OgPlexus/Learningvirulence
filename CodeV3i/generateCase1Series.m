function [Mu, Beta] = generateCase1Series(id)
% generateCase1Series  Deterministic benchmark records of Case 1, with the burn-in removed.
%
% CALLS     dwyerBeta
% CALLED BY runCase1a, runCase1b, runCase1c, buildSensitivitySystems, runSampleSize
%
% INPUT
%   id     one of '1a', '1bi', '1bii', '1biii', '1c'
% OUTPUT
%   Mu, Beta   column vectors of virulence and transmission
switch lower(id)
    case '1a'
        p=2; c_h=0.5; v_star=0.8; b_rho=10; beta_max=1.5; c_g=0.5; k2b=1; a_k=0.5; b_k=3;
        N=6000; burn=100; Mu=zeros(N,1); K2=zeros(N,1); Beta=zeros(N,1);
        Mu(1)=0.30; K2(1)=k2b;
        for n=1:N
            vmax = c_h*(Mu(n)/max(1-Mu(n),1e-14))^(1/p);
            A = exp(1)*vmax/K2(n); Beta(n) = beta_max*A/(A+c_g);
            if n<N
                rho = v_star*b_rho*Mu(n)*exp(1-b_rho*Mu(n));
                Mu(n+1) = rho^p/(rho^p + c_h^p);
                K2(n+1) = k2b*exp(a_k*Beta(n) - b_k*Mu(n));
            end
        end
        Mu=Mu(burn+1:end); Beta=Beta(burn+1:end);

    case '1bi'
        p=2; c_h=0.40; beta_max=1.5; c_g=0.35; t0=1.2; k_v=1.5;
        a=2; b=4; c=15; N=6000; burn=1000; Ntot=N+burn;
        xg = linspace(1e-5,1-1e-5,10000);
        hg = (xg.^a)./(xg.^a + c*(1-xg).^b);
        y = zeros(Ntot,1); y(1)=0.401;
        for n=1:Ntot-1, y(n+1)=4*y(n)*(1-y(n)); end
        V = interp1(hg,xg,y,'pchip'); V = min(max(V,1e-6),1-1e-6);
        Ts = t0./(1+k_v.*V);
        Mu = V.^p./(V.^p + c_h^p);
        A = exp(1).*V.*Ts; Beta = max(beta_max.*A./(A+c_g),0);
        Mu=Mu(burn+1:end); Beta=Beta(burn+1:end);

    case '1bii'
        p=2; c_h=0.40; beta_max=1.5; c_g=0.35; t0=1.2; k_v=1.5;
        N=6000; burn=1000; Ntot=N+burn; ah=1.4; bh=0.3;
        X=zeros(Ntot,1); Y=zeros(Ntot,1); X(1)=0.10; Y(1)=0.30;
        for n=1:Ntot-1, X(n+1)=1-ah*X(n)^2+Y(n); Y(n+1)=bh*X(n); end
        V = (X-min(X))./(max(X)-min(X)); V = 0.05+0.90*V;
        Ts = t0./(1+k_v.*V);
        Mu = V.^p./(V.^p + c_h^p);
        A = exp(1).*V.*Ts; Beta = max(beta_max.*A./(A+c_g),0);
        Mu=Mu(burn+1:end); Beta=Beta(burn+1:end);

    case '1biii'
        beta_max=1.5; N=6000; burn=1000; Ntot=N+burn; ah=1.4; bh=0.3;
        v=zeros(Ntot,1); w=zeros(Ntot,1); v(1)=0.10; w(1)=0.10;
        for n=1:Ntot-1, v(n+1)=1-ah*v(n)^2+w(n); w(n+1)=bh*v(n); end
        v=v(burn+1:end); w=w(burn+1:end);
        Mu   = (v-min(v))/(max(v)-min(v));
        Beta = beta_max*(w-min(w))/(max(w)-min(w));

    case '1c'
        a=2; b=4; c=10; q0=-1.48; q1=0.293; delta=0.01; N=8000; burn=800; Ntot=N+burn;
        xg = linspace(1e-4,1-1e-4,20000)';
        hg = (xg.^a)./(xg.^a + c.*(1-xg).^b);
        Y = zeros(Ntot,1); Y(1)=0.401;
        for n=1:Ntot-1, Y(n+1)=4.*Y(n).*(1-Y(n)); end
        Mu = interp1(hg,xg,Y,'pchip'); Mu = min(max(Mu,0),1);
        Beta = zeros(Ntot,1);
        for n=1:Ntot, Beta(n)=dwyerBeta(Mu(n),delta,q0,q1); end
        Beta = max(Beta,0);
        Mu=Mu(burn+1:end); Beta=Beta(burn+1:end);
end
end
