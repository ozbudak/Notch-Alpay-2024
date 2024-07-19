function [t,mh1Matrix,mh7Matrix]=dCmodelPF(parameters,pb1,pb2)
% Function to simulate the Positive Feedback model with given parameters and perturbations
% Inputs:
%   parameters - Model parameters
%   pb1, pb2 - Perturbation factors for the two cells
% Outputs:
%   t - Time vector
%   mh1Matrix - Matrix of her1 mRNA levels over time
%   mh7Matrix - Matrix of her7 mRNA levels over time

%% Time-Space Settings
L = 2;  % Number of cells
T = 300;  % Total time steps
t = 1:T;  % Time series
%% Model Simulation
% Initialize Matrices to Store Results
mdCMatrix = zeros(L, T);
mdDMatrix = zeros(L, T);
pdCMatrix = zeros(L, T);
pdDMatrix = zeros(L, T);
pdCCMatrix = zeros(L, T);
pdCDMatrix = zeros(L, T);
pdDDMatrix = zeros(L, T);
mh1Matrix = zeros(L, T);
ph1Matrix = zeros(L, T);
mh7Matrix = zeros(L, T);
ph7Matrix = zeros(L, T);
mh6Matrix = zeros(L, T);
ph6Matrix = zeros(L, T);
ph11Matrix = zeros(L, T);
ph77Matrix = zeros(L, T);
ph66Matrix = zeros(L, T);
ph17Matrix = zeros(L, T);
ph16Matrix = zeros(L, T);
ph67Matrix = zeros(L, T);
% Adjust parameters for each cell with perturbations
paraset=zeros(2,64);
paraset(1,:)=parameters.*pb1;
paraset(2,:)=parameters.*pb2;
% Euler's Method Settings
N=100;
h=1/N;
% Loop over each time step
for ti=1:T-1
    % Loop over each cell
    for k=1:L
        if k==1
            parameters=paraset(1,:);
        else
            parameters=paraset(2,:);
        end

        % Model Parameters
        maxmdC_syn = parameters(1);  % Maximum deltaC mRNA synthesis rate
        mdC_deg = parameters(2);  % deltaC mRNA degradation rate
        mdD_syn = parameters(3);  % deltaD mRNA synthesis rate
        mdD_deg = parameters(4);  % deltaD mRNA degradation rate
        pdC_syn = parameters(5);  % DeltaC protein synthesis rate
        pdC_deg = parameters(6);  % DeltaC protein degradation rate
        pdD_syn = parameters(7);  % DeltaD protein synthesis rate
        pdD_deg = parameters(8);  % DeltaD protein degradation rate
        pdCC_deg = parameters(9);  % DeltaC-DeltaC degradation rate
        pdCD_deg = parameters(10);  % DeltaC-DeltaD degradation rate
        pdDD_deg = parameters(11);  % DeltaD-DeltaD degradation rate
        maxmh1_syn = parameters(12);  % Maximum her1 mRNA synthesis rate
        mh1_deg = parameters(13);  % her1 mRNA degradation rate
        ph1_syn = parameters(14);  % Her1 protein synthesis rate
        ph1_deg = parameters(15);  % Her1 protein degradation rate
        maxmh7_syn = parameters(16);  % Maximum her7 mRNA synthesis rate
        mh7_deg = parameters(17);  % her7 mRNA degradation rate
        ph7_syn = parameters(18);  % Her7 protein synthesis rate
        ph7_deg = parameters(19);  % Her7 protein degradation rate
        mh6_syn = parameters(20);  % hes6 mRNA synthesis rate
        mh6_deg = parameters(21);  % hes6 mRNA degradation rate
        ph6_syn = parameters(22);  % Hes6 protein synthesis rate
        ph6_deg = parameters(23);  % Hes6 protein degradation rate
        ph11_deg = parameters(24);  % Her1-Her1 degradation rate
        ph77_deg = parameters(25);  % Her7-Her7 degradation rate
        ph66_deg = parameters(26);  % Hes6-Hes6 degradation rate
        ph17_deg = parameters(27);  % Her1-Her7 degradation rate
        ph16_deg = parameters(28);  % Her1-Hes6 degradation rate
        ph67_deg = parameters(29);  % Hes6-Her7 degradation rate

        pdCC_d = parameters(30);  % Dissociation rate of DeltaC-DeltaC
        pdCC_a = parameters(31);  % Association rate of DeltaC-DeltaC
        pdCD_d = parameters(32);  % Dissociation rate of DeltaC-DeltaD
        pdCD_a = parameters(33);  % Association rate of DeltaC-DeltaD
        pdDD_d = parameters(34);  % Dissociation rate of DeltaD-DeltaD
        pdDD_a = parameters(35);  % Association rate of DeltaD-DeltaD
        ph11_d = parameters(36);  % Dissociation rate of Her1-Her1
        ph11_a = parameters(37);  % Association rate of Her1-Her1
        ph17_d = parameters(38);  % Dissociation rate of Her1-Her7
        ph17_a = parameters(39);  % Association rate of Her1-Her7
        ph16_d = parameters(40);  % Dissociation rate of Her1-Hes6
        ph16_a = parameters(41);  % Association rate of Her1-Hes6
        ph77_d = parameters(42);  % Dissociation rate of Her7-Her7
        ph77_a = parameters(43);  % Association rate of Her7-Her7
        ph67_d = parameters(44);  % Dissociation rate of Hes6-Her7
        ph67_a = parameters(45);  % Association rate of Hes6-Her7
        ph66_d = parameters(46);  % Dissociation rate of Hes6-Hes6
        ph66_a = parameters(47);  % Association rate of Hes6-Hes6

        critph11 = parameters(48);  % DNA-binding dissociation constant for Her1-Her1
        critph67 = parameters(49);  % DNA-binding dissociation constant for Hes6-Her7
        critpdCC = parameters(50);  % DNA-binding dissociation constant for DeltaC-DeltaC
        critpdCD = parameters(51);  % DNA-binding dissociation constant for DeltaC-DeltaD
        critpdDD = parameters(52);  % DNA-binding dissociation constant for DeltaD-DeltaD

        tmdC = round(parameters(53));  % deltaC mRNA transcription time-delay
        tpdC = round(parameters(54));  % DeltaC translation time-delay
        tpdD = round(parameters(55));  % DeltaD translation time-delay
        tmh1 = round(parameters(56));  % her1 mRNA transcription time-delay
        tph1 = round(parameters(57));  % Her1 translation time-delay
        tmh7 = round(parameters(58));  % her7 mRNA transcription time-delay
        tph7 = round(parameters(59));  % Her7 translation time-delay
        tph6 = round(parameters(60));  % Hes6 translation time-delay

        E1 = 1 / (1 + exp(parameters(61)));  % Efficiency parameter of constant transcription
        ECC = 1 / (1 + exp(parameters(61) - parameters(62)));  % Efficiency parameter of transcription activated by DeltaC-DeltaC
        ECD = 1 / (1 + exp(parameters(61) - parameters(63)));  % Efficiency parameter of transcription activated by DeltaC-DeltaD
        EDD = 1 / (1 + exp(parameters(61) - parameters(64)));  % Efficiency parameter of transcription activated by DeltaD-DeltaD

        mdC=zeros(N+1,1);
        mdD=zeros(N+1,1);
        pdC=zeros(N+1,1);
        pdD=zeros(N+1,1);
        pdCC=zeros(N+1,1);
        pdCD=zeros(N+1,1);
        pdDD=zeros(N+1,1);
        mh1=zeros(N+1,1);
        ph1=zeros(N+1,1);
        mh7=zeros(N+1,1);
        ph7=zeros(N+1,1);
        mh6=zeros(N+1,1);
        ph6=zeros(N+1,1);
        ph11=zeros(N+1,1);
        ph77=zeros(N+1,1);
        ph66=zeros(N+1,1);
        ph17=zeros(N+1,1);
        ph16=zeros(N+1,1);
        ph67=zeros(N+1,1);

        mdC(1)=mdCMatrix(k,ti);
        mdD(1)=mdDMatrix(k,ti);
        pdC(1)=pdCMatrix(k,ti);
        pdD(1)=pdDMatrix(k,ti);
        pdCC(1)=pdCCMatrix(k,ti);
        pdCD(1)=pdCDMatrix(k,ti);
        pdDD(1)=pdDDMatrix(k,ti);
        mh1(1)=mh1Matrix(k,ti);
        ph1(1)=ph1Matrix(k,ti);
        mh7(1)=mh7Matrix(k,ti);
        ph7(1)=ph7Matrix(k,ti);
        mh6(1)=mh6Matrix(k,ti);
        ph6(1)=ph6Matrix(k,ti);
        ph11(1)=ph11Matrix(k,ti);
        ph77(1)=ph77Matrix(k,ti);
        ph66(1)=ph66Matrix(k,ti);
        ph17(1)=ph17Matrix(k,ti);
        ph16(1)=ph16Matrix(k,ti);
        ph67(1)=ph67Matrix(k,ti);

        % Euler's method
        for n=1:N
            % Update deltaC mRNA (With Positive Feedback)
            if ti>tmdC
                if k==1
                    den=1+(ph11Matrix(k,ti-tmdC)/critph11)^2+(ph67Matrix(k,ti-tmdC)/critph67)^2+(pdCCMatrix(2,ti-tmdC)/critpdCC)+(pdCDMatrix(2,ti-tmdC)/critpdCD)+(pdDDMatrix(2,ti-tmdC)/critpdDD);
                    mdC(n+1)=mdC(n)+h*(maxmdC_syn*(1/den*E1+(pdCCMatrix(2,ti-tmdC)/critpdCC)/den*ECC+(pdCDMatrix(2,ti-tmdC)/critpdCD)/den*ECD+(pdDDMatrix(2,ti-tmdC)/critpdDD)/den*EDD)-mdC_deg*mdC(n));
                else
                    den=1+(ph11Matrix(k,ti-tmdC)/critph11)^2+(ph67Matrix(k,ti-tmdC)/critph67)^2+(pdCCMatrix(1,ti-tmdC)/critpdCC)+(pdCDMatrix(1,ti-tmdC)/critpdCD)+(pdDDMatrix(1,ti-tmdC)/critpdDD);
                    mdC(n+1)=mdC(n)+h*(maxmdC_syn*(1/den*E1+(pdCCMatrix(1,ti-tmdC)/critpdCC)/den*ECC+(pdCDMatrix(1,ti-tmdC)/critpdCD)/den*ECD+(pdDDMatrix(1,ti-tmdC)/critpdDD)/den*EDD)-mdC_deg*mdC(n));
                end
            else
                mdC(n+1)=mdC(n);
            end

            % Update deltaD mRNA
            mdD(n+1)=mdD(n)+h*(mdD_syn-mdD_deg*mdD(n));

            % Update DeltaC protein
            if ti>tpdC
                pdC(n+1)=pdC(n)+h*(pdC_syn*mdCMatrix(k,ti-tpdC)-pdC_deg*pdC(n)+2*pdCC_d*pdCC(n)+pdCD_d*pdCD(n)-2*pdCC_a*pdC(n)*pdC(n)-pdCD_a*pdC(n)*pdD(n));
            else
                pdC(n+1)=pdC(n);
            end

            % Update DeltaD protein
            if ti>tpdD
                pdD(n+1)=pdD(n)+h*(pdD_syn*mdDMatrix(k,ti-tpdD)-pdD_deg*pdD(n)+2*pdDD_d*pdDD(n)+pdCD_d*pdCD(n)-2*pdDD_a*pdD(n)*pdD(n)-pdCD_a*pdC(n)*pdD(n));
            else
                pdD(n+1)=pdD(n);
            end

            % Update DeltaC-DeltaC dimer
            pdCC(n+1)=pdCC(n)+h*(pdCC_a*pdC(n)*pdC(n)-pdCC_d*pdCC(n)-pdCC_deg*pdCC(n));
            % Update DeltaC-DeltaD dimer
            pdCD(n+1)=pdCD(n)+h*(pdCD_a*pdC(n)*pdD(n)-pdCD_d*pdCD(n)-pdCD_deg*pdCD(n));
            % Update DeltaD-DeltaD dimer
            pdDD(n+1)=pdDD(n)+h*(pdDD_a*pdD(n)*pdD(n)-pdDD_d*pdDD(n)-pdDD_deg*pdDD(n));

            % Update her1 mRNA
            if ti>tmh1
                if k==1
                    den=1+(ph11Matrix(k,ti-tmh1)/critph11)^2+(ph67Matrix(k,ti-tmh1)/critph67)^2+(pdCCMatrix(2,ti-tmh1)/critpdCC)+(pdCDMatrix(2,ti-tmh1)/critpdCD)+(pdDDMatrix(2,ti-tmh1)/critpdDD);
                    mh1(n+1)=mh1(n)+h*(maxmh1_syn*(1/den*E1+(pdCCMatrix(2,ti-tmh1)/critpdCC)/den*ECC+(pdCDMatrix(2,ti-tmh1)/critpdCD)/den*ECD+(pdDDMatrix(2,ti-tmh1)/critpdDD)/den*EDD)-mh1_deg*mh1(n));
                else
                    den=1+(ph11Matrix(k,ti-tmh1)/critph11)^2+(ph67Matrix(k,ti-tmh1)/critph67)^2+(pdCCMatrix(1,ti-tmh1)/critpdCC)+(pdCDMatrix(1,ti-tmh1)/critpdCD)+(pdDDMatrix(1,ti-tmh1)/critpdDD);
                    mh1(n+1)=mh1(n)+h*(maxmh1_syn*(1/den*E1+(pdCCMatrix(1,ti-tmh1)/critpdCC)/den*ECC+(pdCDMatrix(1,ti-tmh1)/critpdCD)/den*ECD+(pdDDMatrix(1,ti-tmh1)/critpdDD)/den*EDD)-mh1_deg*mh1(n));
                end
            else
                mh1(n+1)=mh1(n);
            end

            % Update Her1 protein
            if ti>tph1
                ph1(n+1)=ph1(n)+h*(ph1_syn*mh1Matrix(k,ti-tph1)-ph1_deg*ph1(n)+2*ph11_d*ph11(n)+ph17_d*ph17(n)+ph16_d*ph16(n)-2*ph11_a*ph1(n)*ph1(n)-ph17_a*ph1(n)*ph7(n)-ph16_a*ph1(n)*ph6(n));
            else
                ph1(n+1)=ph1(n);
            end

            % Update her7 mRNA
            if ti>tmh7
                if k==1
                    den=1+(ph11Matrix(k,ti-tmh7)/critph11)^2+(ph67Matrix(k,ti-tmh7)/critph67)^2+(pdCCMatrix(2,ti-tmh7)/critpdCC)+(pdCDMatrix(2,ti-tmh7)/critpdCD)+(pdDDMatrix(2,ti-tmh7)/critpdDD);
                    mh7(n+1)=mh7(n)+h*(maxmh7_syn*(1/den*E1+(pdCCMatrix(2,ti-tmh7)/critpdCC)/den*ECC+(pdCDMatrix(2,ti-tmh7)/critpdCD)/den*ECD+(pdDDMatrix(2,ti-tmh7)/critpdDD)/den*EDD)-mh7_deg*mh7(n));
                else
                    den=1+(ph11Matrix(k,ti-tmh7)/critph11)^2+(ph67Matrix(k,ti-tmh7)/critph67)^2+(pdCCMatrix(1,ti-tmh7)/critpdCC)+(pdCDMatrix(1,ti-tmh7)/critpdCD)+(pdDDMatrix(1,ti-tmh7)/critpdDD);
                    mh7(n+1)=mh7(n)+h*(maxmh7_syn*(1/den*E1+(pdCCMatrix(1,ti-tmh7)/critpdCC)/den*ECC+(pdCDMatrix(1,ti-tmh7)/critpdCD)/den*ECD+(pdDDMatrix(1,ti-tmh7)/critpdDD)/den*EDD)-mh7_deg*mh7(n));
                end
            else
                mh7(n+1)=mh7(n);
            end

            % Update Her7 protein
            if ti>tph7
                ph7(n+1)=ph7(n)+h*(ph7_syn*mh7Matrix(k,ti-tph7)-ph7_deg*ph7(n)+2*ph77_d*ph77(n)+ph17_d*ph17(n)+ph67_d*ph67(n)-2*ph77_a*ph7(n)*ph7(n)-ph17_a*ph1(n)*ph7(n)-ph67_a*ph6(n)*ph7(n));
            else
                ph7(n+1)=ph7(n);
            end

            % Update hes6 mRNA
            mh6(n+1)=mh6(n)+h*(mh6_syn-mh6_deg*mh6(n));

            % Update Hes6 protein
            if ti>tph6
                ph6(n+1)=ph6(n)+h*(ph6_syn*mh6Matrix(k,ti-tph6)-ph6_deg*ph6(n)+2*ph66_d*ph66(n)+ph16_d*ph16(n)+ph67_d*ph67(n)-2*ph66_a*ph6(n)*ph6(n)-ph16_a*ph1(n)*ph6(n)-ph67_a*ph6(n)*ph7(n));
            else
                ph6(n+1)=ph6(n);
            end

            % Update Her1-Her1 dimer
            ph11(n+1)=ph11(n)+h*(ph11_a*ph1(n)*ph1(n)-ph11_d*ph11(n)-ph11_deg*ph11(n));
            % Update Her7-Her7 dimer
            ph77(n+1)=ph77(n)+h*(ph77_a*ph7(n)*ph7(n)-ph77_d*ph77(n)-ph77_deg*ph77(n));
            % Update Hes6-Hes6 dimer
            ph66(n+1)=ph66(n)+h*(ph66_a*ph6(n)*ph6(n)-ph66_d*ph66(n)-ph66_deg*ph66(n));
            % Update Her1-Her7 dimer
            ph17(n+1)=ph17(n)+h*(ph17_a*ph1(n)*ph7(n)-ph17_d*ph17(n)-ph17_deg*ph17(n));
            % Update Her1-Hes6 dimer
            ph16(n+1)=ph16(n)+h*(ph16_a*ph1(n)*ph6(n)-ph16_d*ph16(n)-ph16_deg*ph16(n));
            % Update Hes-Her7 dimer
            ph67(n+1)=ph67(n)+h*(ph67_a*ph6(n)*ph7(n)-ph67_d*ph67(n)-ph67_deg*ph67(n));
        end
        
        % Store the results for the next time step
        mdCMatrix(k,ti+1)=mdC(N+1);
        mdDMatrix(k,ti+1)=mdD(N+1);
        pdCMatrix(k,ti+1)=pdC(N+1);
        pdDMatrix(k,ti+1)=pdD(N+1);
        pdCCMatrix(k,ti+1)=pdCC(N+1);
        pdCDMatrix(k,ti+1)=pdCD(N+1);
        pdDDMatrix(k,ti+1)=pdDD(N+1);
        mh1Matrix(k,ti+1)=mh1(N+1);
        ph1Matrix(k,ti+1)=ph1(N+1);
        mh7Matrix(k,ti+1)=mh7(N+1);
        ph7Matrix(k,ti+1)=ph7(N+1);
        mh6Matrix(k,ti+1)=mh6(N+1);
        ph6Matrix(k,ti+1)=ph6(N+1);
        ph11Matrix(k,ti+1)=ph11(N+1);
        ph77Matrix(k,ti+1)=ph77(N+1);
        ph66Matrix(k,ti+1)=ph66(N+1);
        ph17Matrix(k,ti+1)=ph17(N+1);
        ph16Matrix(k,ti+1)=ph16(N+1);
        ph67Matrix(k,ti+1)=ph67(N+1);
    end
end
end