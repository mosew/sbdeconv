function y = delaymodel(a0,ell,u,t)
    global N n
    MN = build_MN(ell);
    n = length(t);
    y = zeros(1,n);
    assert(size(t,1)==1);

    KN = build_KN(a0);
    AN = MN\KN;
    BN = [zeros(1,N+1),2*N/ell,zeros(1,N-1)]';
    CN = [1,zeros(1,2*N)]*MN;
    
    if length(u)==1
        model = @(r,u) integral(@(s) model0(r-s,AN,BN,CN).*u(s), 0,r);
        for i = 1:n
            y(i) = model(t(i),u);
        end
    else
        for j = 2:n
            for k = 1:j
                y(j) = y(j)+CN*sg(t(j)-t(k),AN)*BN*u(k);
            end
        end
    end
end

function gg = sg(s,AN)
    global N
    gg = zeros(2*N+1,2*N+1,length(s));
    for i = 1:length(s)
        gg(:,:,i) = expm(s(i)*AN);
    end
end

function m = model0(s,AN,BN,CN)
    m=zeros(size(s));
    for i = 1:length(s)
        m(i) = CN*sg(s(i),AN)*BN;
    end
end