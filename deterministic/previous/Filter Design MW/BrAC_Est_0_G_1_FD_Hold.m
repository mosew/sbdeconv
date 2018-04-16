function [Est_BrAC_TAC] = BrAC_Est_0_G_1_FD(TAC,h_In)


%  Input:  2-Column Array TAC:  First Column: Time in Hours
%                               Second Column: TAC in % Alcohol
%          Row Vector h_In: h_In(1): r_1
%                           h_In(2): r_2
%                           h_In(3:end): Convolution Kernel in units of % 
%                                        Alcohol vs time in hours sampled 
%                                        at 1 minute (1/60 hours)time 
%                                        intervals

%  Output:  2-Column Array Est_BrAC:  First Column: Time in Hours
%                                     Second Column: Estimated BrAC 
%                                     in % Alcohol

% r_1_Default, r_2_Default, and h_Default were computed using the
% Brown_Data_lab data using R21_BrAC_Estimator_L_Paper_F_MM_NB_A_Brown_C
% and using all seven cases (1:7).

n_index = 6;


r_1_Default = 4.7733;	
r_2_Default = 1.7020;

h_Default = [1.08E-16	0	1.26E-16	4.83E-15	0	0	1.87E-09	1.92E-07	2.51E-06	1.51E-05	5.81E-05 ...
    0.00016685	0.00039036	0.00078616	0.0014132	0.0023256	0.0035674	0.0051702	0.0071523	0.0095188...
    0.012263	0.01537	0.018816	0.022571	0.026603	0.030877	0.035357	0.040009	0.044797	0.049689...
    0.054655	0.059665	0.064695	0.069719	0.074718	0.079673	0.084566	0.089385	0.094115	0.098748...
    0.10327	0.10768	0.11198	0.11614	0.12018	0.12409	0.12786	0.1315	0.13501	0.13838	0.14162	0.14473	0.14771	0.15057...
    0.15329	0.1559	0.15838	0.16074	0.163	0.16513	0.16716	0.16909	0.17091	0.17263	0.17426	0.17579	0.17724	0.17859...
    0.17986	0.18106	0.18217	0.1832	0.18417	0.18506	0.18589	0.18665	0.18734	0.18798	0.18856	0.18908	0.18955	0.18997...
    0.19034	0.19066	0.19094	0.19117	0.19136	0.19151	0.19163	0.1917	0.19174	0.19175	0.19172	0.19166	0.19157	0.19146...
    0.19131	0.19114	0.19095	0.19073	0.19049	0.19023	0.18994	0.18964	0.18932	0.18898	0.18862	0.18825	0.18786	0.18745...
    0.18703	0.1866	0.18615	0.1857	0.18523	0.18475	0.18425	0.18375	0.18324	0.18272	0.18219	0.18166	0.18111	0.18056...
    0.18001	0.17944	0.17887	0.1783	0.17771	0.17713	0.17654	0.17594	0.17534	0.17474	0.17413	0.17352	0.17291	0.1723...
    0.17168	0.17106	0.17043	0.16981	0.16918	0.16855	0.16792	0.16729	0.16666	0.16603	0.1654	0.16476	0.16413	0.16349...
    0.16286	0.16222	0.16159	0.16095	0.16032	0.15968	0.15905	0.15841	0.15778	0.15715	0.15651	0.15588	0.15525	0.15462...
    0.15399	0.15337	0.15274	0.15212	0.15149	0.15087	0.15025	0.14963	0.14901	0.14839	0.14778	0.14716	0.14655	0.14594...
    0.14533	0.14473	0.14412	0.14352	0.14292	0.14232	0.14172	0.14112	0.14053	0.13993	0.13934	0.13875	0.13817	0.13758...
    0.137	0.13642	0.13584	0.13526	0.13469	0.13412	0.13354	0.13298	0.13241	0.13184	0.13128	0.13072	0.13016	0.12961...
    0.12905	0.1285	0.12795	0.1274	0.12686	0.12631	0.12577	0.12523	0.1247	0.12416	0.12363	0.1231	0.12257	0.12204...
    0.12152	0.12099	0.12047	0.11995	0.11944	0.11892	0.11841	0.1179	0.11739	0.11689	0.11638	0.11588	0.11538	0.11488...
    0.11439	0.11389	0.1134	0.11291	0.11243	0.11194	0.11146	0.11098	0.1105	0.11002	0.10954	0.10907	0.1086	0.10813...
    0.10766	0.1072	0.10673	0.10627	0.10581	0.10535	0.1049	0.10444	0.10399	0.10354	0.10309	0.10265	0.1022	0.10176...
    0.10132	0.10088	0.10045	0.10001	0.099578	0.099147	0.098718	0.09829	0.097865	0.097441	0.097019...
    0.096599	0.09618	0.095764	0.095349	0.094936	0.094525	0.094115	0.093707	0.093302	0.092897...
    0.092495	0.092094	0.091695	0.091298	0.090902	0.090508	0.090116	0.089725	0.089337	0.088949...
    0.088564	0.08818	0.087798	0.087417	0.087038	0.086661	0.086286	0.085912	0.085539	0.085168...
    0.084799	0.084432	0.084066	0.083701	0.083338	0.082977	0.082617	0.082259	0.081903	0.081547...
    0.081194	0.080842	0.080491	0.080143	0.079795	0.079449	0.079105	0.078762	0.07842	0.07808	0.077742...
    0.077405	0.077069	0.076735	0.076402	0.076071	0.075741	0.075413	0.075086	0.07476	0.074436...
    0.074113	0.073792	0.073472	0.073153	0.072836	0.07252	0.072206	0.071893	0.071581	0.071271...
    0.070961	0.070654	0.070347	0.070042	0.069739	0.069436	0.069135	0.068835	0.068537	0.06824...
    0.067944	0.067649	0.067356	0.067064	0.066773	0.066483	0.066195	0.065908	0.065622	0.065338...
    0.065054	0.064772	0.064491	0.064212	0.063933	0.063656	0.06338	0.063105	0.062831	0.062559...
    0.062288	0.062017	0.061749	0.061481	0.061214	0.060949	0.060684	0.060421	0.060159	0.059898...
    0.059639	0.05938	0.059122	0.058866	0.058611	0.058357	0.058104	0.057852	0.057601	0.057351...
    0.057102	0.056855	0.056608	0.056363	0.056118	0.055875	0.055632	0.055391	0.055151	0.054912...
    0.054674	0.054437	0.054201	0.053965	0.053731	0.053498	0.053266	0.053035	0.052805	0.052576...
    0.052348	0.052121	0.051895	0.05167	0.051446	0.051223	0.051001	0.05078	0.05056	0.05034	0.050122...
    0.049905	0.049688	0.049473	0.049258	0.049045	0.048832	0.04862	0.048409	0.048199	0.04799...
    0.047782	0.047575	0.047369	0.047163	0.046959	0.046755	0.046552	0.046351	0.04615	0.045949...
    0.04575	0.045552	0.045354	0.045157	0.044962	0.044767	0.044573	0.044379	0.044187	0.043995...
    0.043804	0.043614	0.043425	0.043237	0.043049	0.042863	0.042677	0.042492	0.042308	0.042124...
    0.041941	0.041759	0.041578	0.041398	0.041219	0.04104	0.040862	0.040685	0.040508	0.040333...
    0.040158	0.039983	0.03981	0.039637	0.039466	0.039294	0.039124	0.038954	0.038785	0.038617...
    0.03845	0.038283	0.038117	0.037952	0.037787	0.037623	0.03746	0.037298	0.037136	0.036975...
    0.036814	0.036655	0.036496	0.036338	0.03618	0.036023	0.035867	0.035711	0.035556	0.035402...
    0.035249	0.035096	0.034944	0.034792	0.034641	0.034491	0.034341	0.034193	0.034044	0.033897...
    0.03375	0.033603	0.033458	0.033312	0.033168	0.033024	0.032881	0.032738	0.032596	0.032455...
    0.032314	0.032174	0.032035	0.031896	0.031757	0.03162	0.031482	0.031346	0.03121	0.031075...	
    0.03094	0.030806	0.030672	0.030539	0.030407	0.030275	0.030144	0.030013	0.029883	0.029753...
    0.029624	0.029496	0.029368	0.02924	0.029113	0.028987	0.028862	0.028736	0.028612	0.028488...
    0.028364	0.028241	0.028119	0.027997	0.027875	0.027754	0.027634	0.027514	0.027395	0.027276...
    0.027158	0.02704	0.026923	0.026806	0.02669	0.026574	0.026459	0.026344	0.02623	0.026116	0.026003...
    0.02589	0.025778	0.025666	0.025555	0.025444	0.025333	0.025224	0.025114	0.025005	0.024897...
    0.024789	0.024681	0.024574	0.024468	0.024362	0.024256	0.024151	0.024046	0.023942	0.023838...
    0.023735	0.023632	0.023529	0.023427	0.023326	0.023224	0.023124	0.023023	0.022924	0.022824...
    0.022725	0.022627	0.022529	0.022431	0.022334	0.022237	0.02214	0.022044	0.021949	0.021853...
    0.021759	0.021664	0.02157	0.021477	0.021384	0.021291	0.021199	0.021107	0.021015	0.020924...
    0.020833	0.020743	0.020653	0.020563	0.020474	0.020385	0.020297	0.020209	0.020121	0.020034...
    0.019947	0.019861	0.019775	0.019689	0.019603	0.019518	0.019434	0.01935	0.019266	0.019182...
    0.019099	0.019016	0.018934	0.018851	0.01877	0.018688	0.018607	0.018527	0.018446	0.018366...
    0.018287	0.018207	0.018128	0.01805	0.017971	0.017894	0.017816	0.017739	0.017662	0.017585...
    0.017509	0.017433	0.017357	0.017282	0.017207	0.017133	0.017058	0.016984	0.016911	0.016837...
    0.016764	0.016692	0.016619	0.016547	0.016475	0.016404	0.016333	0.016262	0.016191	0.016121...
    0.016051	0.015982	0.015912	0.015843	0.015775	0.015706	0.015638	0.01557	0.015503	0.015436...
    0.015369	0.015302	0.015236	0.01517	0.015104	0.015038	0.014973	0.014908	0.014843	0.014779...
    0.014715	0.014651	0.014588	0.014524	0.014461	0.014399	0.014336	0.014274	0.014212	0.014151...
    0.014089	0.014028	0.013967	0.013907	0.013846	0.013786	0.013726	0.013667	0.013608	0.013549...
    0.01349	0.013431	0.013373	0.013315	0.013257	0.0132	0.013143	0.013086	0.013029	0.012972	0.012916...
    0.01286	0.012804	0.012749	0.012694	0.012639	0.012584	0.012529	0.012475	0.012421	0.012367	0.012313...
    0.01226	0.012207	0.012154	0.012101	0.012049	0.011996	0.011944	0.011892	0.011841	0.01179	0.011738...
    0.011688	0.011637	0.011586	0.011536	0.011486	0.011436	0.011387	0.011337	0.011288	0.011239...
    0.01119	0.011142	0.011094	0.011046	0.010998	0.01095	0.010902	0.010855	0.010808	0.010761	0.010715...
    0.010668	0.010622	0.010576	0.01053	0.010484	0.010439	0.010393	0.010348	0.010304	0.010259	0.010214...
    0.01017	0.010126	0.010082	0.010038	0.0099948	0.0099514	0.0099083	0.0098653	0.0098225	0.0097799	0.0097375...
    0.0096953	0.0096533	0.0096114	0.0095697	0.0095282	0.0094869	0.0094457	0.0094048	0.009364	0.0093234	0.009283...
    0.0092427	0.0092026	0.0091627	0.009123	0.0090834	0.009044	0.0090048	0.0089657	0.0089269	0.0088882	0.0088496...
    0.0088112	0.008773	0.008735	0.0086971	0.0086594	0.0086218	0.0085844	0.0085472	0.0085101	0.0084732	0.0084365...
    0.0083999	0.0083635	0.0083272	0.0082911	0.0082551	0.0082193	0.0081837	0.0081482	0.0081129	0.0080777	0.0080427...
    0.0080078	0.007973	0.0079385	0.007904	0.0078698	0.0078356	0.0078017	0.0077678	0.0077341	0.0077006	0.0076672...
    0.007634	0.0076008	0.0075679	0.0075351	0.0075024	0.0074699	0.0074375	0.0074052	0.0073731	0.0073411	0.0073093...
    0.0072776	0.007246	0.0072146	0.0071833	0.0071522	0.0071211	0.0070903	0.0070595	0.0070289	0.0069984	0.0069681...
    0.0069379	0.0069078	0.0068778	0.006848	0.0068183	0.0067887	0.0067593	0.00673	0.0067008	0.0066717	0.0066428...
    0.006614	0.0065853	0.0065567	0.0065283	0.0065	0.0064718	0.0064437	0.0064158	0.006388	0.0063603	0.0063327...
    0.0063052	0.0062779	0.0062507	0.0062235	0.0061966	0.0061697	0.0061429	0.0061163	0.0060898	0.0060634...
    0.0060371	0.0060109	0.0059848	0.0059589	0.005933	0.0059073	0.0058817	0.0058562	0.0058308	0.0058055...
    0.0057803	0.0057552	0.0057303	0.0057054	0.0056807	0.0056561	0.0056315	0.0056071	0.0055828	0.0055586...
    0.0055345	0.0055105	0.0054866	0.0054628	0.0054391	0.0054155	0.005392	0.0053686	0.0053453	0.0053222	0.0052991...
    0.0052761	0.0052532	0.0052304	0.0052078	0.0051852	0.0051627	0.0051403	0.005118	0.0050958	0.0050737	0.0050517...
    0.0050298	0.005008	0.0049863	0.0049647	0.0049431	0.0049217	0.0049003	0.0048791	0.0048579	0.0048369	0.0048159...
    0.004795	0.0047742	0.0047535	0.0047329	0.0047124	0.0046919	0.0046716	0.0046513	0.0046312	0.0046111	0.0045911...
    0.0045712	0.0045513	0.0045316	0.0045119	0.0044924	0.0044729	0.0044535	0.0044342	0.004415	0.0043958	0.0043768...
    0.0043578	0.0043389	0.0043201	0.0043013	0.0042827	0.0042641	0.0042456	0.0042272	0.0042089	0.0041906	0.0041724...
    0.0041543	0.0041363	0.0041184	0.0041005	0.0040827	0.004065	0.0040474	0.0040299	0.0040124	0.003995	0.0039777...
    0.0039604	0.0039432	0.0039261	0.0039091	0.0038922	0.0038753	0.0038585	0.0038417	0.0038251	0.0038085	0.003792...
    0.0037755	0.0037592	0.0037429	0.0037266	0.0037105	0.0036944	0.0036783	0.0036624	0.0036465	0.0036307	0.003615...
    0.0035993	0.0035837	0.0035681	0.0035527	0.0035372	0.0035219	0.0035066	0.0034914	0.0034763	0.0034612	0.0034462...
    0.0034313	0.0034164	0.0034016	0.0033868	0.0033721	0.0033575	0.0033429	0.0033284	0.003314	0.0032996	0.0032853...
    0.0032711	0.0032569	0.0032428	0.0032287	0.0032147	0.0032008	0.0031869	0.0031731	0.0031593	0.0031456	0.003132...
    0.0031184	0.0031049	0.0030914	0.003078	0.0030646	0.0030513	0.0030381	0.0030249	0.0030118	0.0029988	0.0029857...
    0.0029728	0.0029599	0.0029471	0.0029343	0.0029216	0.0029089	0.0028963	0.0028837	0.0028712	0.0028588	0.0028464...
    0.002834	0.0028217	0.0028095	0.0027973	0.0027852	0.0027731	0.0027611	0.0027491	0.0027372	0.0027253	0.0027135...
    0.0027017	0.00269	0.0026783	0.0026667	0.0026552	0.0026436	0.0026322	0.0026208	0.0026094	0.0025981	0.0025868...
    0.0025756	0.0025644	0.0025533	0.0025422	0.0025312	0.0025202	0.0025093	0.0024984	0.0024876	0.0024768	0.0024661...
    0.0024554	0.0024447	0.0024341	0.0024236	0.0024131	0.0024026	0.0023922	0.0023818	0.0023715	0.0023612	0.0023509...
    0.0023407	0.0023306	0.0023205	0.0023104	0.0023004	0.0022904	0.0022805	0.0022706	0.0022608	0.002251	0.0022412...
    0.0022315	0.0022218	0.0022122	0.0022026	0.002193	0.0021835	0.002174	0.0021646	0.0021552	0.0021459	0.0021366...
    0.0021273	0.0021181	0.0021089	0.0020997	0.0020906	0.0020816	0.0020725	0.0020636	0.0020546	0.0020457	0.0020368...
    0.002028	0.0020192	0.0020104	0.0020017	0.001993	0.0019844	0.0019758	0.0019672	0.0019587	0.0019502	0.0019417...
    0.0019333	0.0019249	0.0019166	0.0019083	0.0019	0.0018918	0.0018836	0.0018754	0.0018673	0.0018592	0.0018511...
    0.0018431	0.0018351	0.0018271	0.0018192	0.0018113	0.0018035	0.0017956	0.0017878	0.0017801	0.0017724	0.0017647...
    0.001757	0.0017494	0.0017418	0.0017343];

if (norm(h_In) ~= 0)
    
    r_1 = h_In(1);
    r_2 = h_In(2);
    h = h_In(3:end);
    
else
    
    r_1 = r_1_Default;
    r_2 = r_2_Default;
    h = h_Default;
    
end

splines_per_hour = [1,2,3,4,5,6,10,12,15,20,30,60];
n = splines_per_hour(n_index);
m = 60/n;
phi_n = [(0:m)/m,((m-1):-1:0)/m];

T_h = (length(h)-1)/60;
T_y = TAC(end,1);

N = floor(n*T_y);

if (N*m > 60*T_h)
    
    h = [h,zeros(1,N*m - 60*T_h)];
    
end  

A = zeros(N*m + 1,N);

for i = 1:N*m
    
   j = 1;
   
   while (m*(j-1)<i)
       
       k = fix(j/n);
       l = rem(j,n);
       
        A(i+1,j) = 0;
        
        p_min = max(0,60*k + (l-1)*m);
        p_max = min(i,60*k + (l+1)*m);
   
        for p = p_min:p_max
            
            A(i+1,j) = A(i+1,j) + (1/60)*h(1 + i-p)*phi_n(1 + p - p_min);
            
        end
        
        j = j + 1;
           
   end
   
end

D = [(1/N)*[(2/3)*ones(1,N-1),(1/3)]];  
S = [(1/N)*[(1/6)*ones(1,N-1)]];   
Q_M = (N*m/60)*(diag(D) + diag(S,1) + diag(S,-1));

D = [N*[2*ones(1,N-1),1]]; 
S = [-N*[ones(1,N-1)]];   
Q_K = (60/(N*m))*(diag(D) + diag(S,1) + diag(S,-1));

time = (0:N*m)'/60;
y = max(spline(TAC(:,1),1000*TAC(:,2),time),0);

A_Tilda = [A;r_1*Q_M + r_2*Q_K];
b_Tilda = [y;zeros(N,1)];

OPTIONS = optimset('Display','off');
warning off
eBrAC = lsqnonneg(A_Tilda,b_Tilda,OPTIONS);
warning on

eTAC = A*eBrAC/1000;

Est_BrAC_TAC= [time,max(spline((0:m:N*m)'/60,.001*[0;eBrAC],time),0),eTAC];
end

