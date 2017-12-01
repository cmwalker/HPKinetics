classdef AbstractPools 
    %ABSTRACTPOOLS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function defaults()
            % DEFAULTS explains the default values for each parameter
            names = {'ExchangeTerms','T1s','FaList','TRList',...
                'PerfusionTerms','volumeFractions','VIF','fitOptions'};
            discriptions = {'A  NxN Matrix of Exchange Terms, where N is the number of chemical pools. The From pools should be along the rRows With the To pool along the Columns. Diagnal elemets will be set to zero'...
                ' A  Row vector of T1 decay times for each chemical pool.'...
                ' A  NxM of matrix of flip angles in radians, where N is the number of excitations and M is the number of chemical Pools'...
                ' A  NxM of Excitation Times in seconds, where N is the number of excitations and M is the number of chemical Pools'...
                ' A Row Vector of perfusion Exchange Constnats for each chemical pool.'...
                ' A Row Vector of volme fraction for each chemical pool. Only one value can be use if all pools have the same volume fraction.'...
                ' A function of a time variable (t) in seconds that returns a Row vector for the VIF of each chemical pool at the time t.'...
                ' Matlab FitOptions object'};
            defaultsVals = {'0','100','0','0','0','1','@(t)0','optimset(''lsqcurvefit'')'};
            fprintf('*Note* all terms must be a vector of size 1 x N where N is the number of chemical Pools\n')
            for i = 1:numel(names)
                fprintf('''%s'': %s\n Default Vaule: %s\n',...
                    names{i},discriptions{i},defaultsVals{i});
            end
        end
        function paramsOut = parseParams(paramsIn)
            % parseParams: a function to fill default param values if they are
            % not defined
            default = struct('ExchangeTerms',0,'T1s',100,'FaList',0,...
                'TRList',0,'PerfusionTerms',0,'volumeFractions',1,...
                'fitOptions', optimset('lsqcurvefit'));
            tmpNames = fieldnames(default);
            paramsOut = paramsIn;
            for i = 1:numel(tmpNames)
                if ~isfield(paramsOut,tmpNames{i})
                    paramsOut.(tmpNames{i}) = default.(tmpNames{i});
                end
            end
            P = paramsOut.PerfusionTerms; % Matrix of Physical Exchange Terms
            C = paramsOut.ExchangeTerms; % Matrix of Chemical Exchange Terms
            V = paramsOut.volumeFractions; % Vector of Volume Fractions
            T1 = paramsIn.T1s;
            nC = size(P,1);
            nP = size(C,1);
            % Fill all flip angles with a value if only one flip angle is passed in
            if size(paramsOut.FaList,2)==1
                paramsOut.FaList = repmat(paramsOut.FaList(:,1),...
                    1,length(paramsOut.TRList));
            end
            paramsOut.FaList = repmat(paramsOut.FaList,nC*nP/size(paramsOut.FaList,1),1);
            % Fill the volume fractions
            tmpCase = nP - numel(V);
            switch tmpCase
                case 1
                    V(end+1) = 1-...
                        sum(V);
                    if sum(V) > 1
                        error('volume fractions added up to more than 1.\n')
                    end
                case 0
%                     V = repelem(...
%                         V,nC);
                otherwise
                    error('To many or not enough volume fractions passed in.\n')
            end
            paramsOut.volumeFractions = V;
            % Create T1 decay matrix
            if numel(T1) == nC
                if size(T1,1) ~= nC
                    T1 = T1.';
                end
                T1 = repmat(T1,[nP 1]);
            end
            % Create A matrix
            tmpA = zeros(nP,nP,nC,nC);
            A = zeros(nP*nC,nP*nC);
            R1 = diag(1./T1);
            for i = 1:nP
                for j = 1:nP
                    if i == j % Set up Diagnal Elements
                        tmpA(i,j,:,:) = squeeze(C(i,:,:))-...
                            diag(sum(squeeze(C(i,:,:)),1)+...
                            sum(squeeze(P(:,:,i)),2).'./V(i));
                    else % Set up off diagonal elements
                        tmpA(i,j,:,:) = diag(P(:,i,j))./V(i);
                    end
                    % Reduce into final A matirx
                    A(i*nC-nC+1:i*nC,j*nC-nC+1:j*nC) = tmpA(i,j,:,:);
                end
            end
            A = A-R1; % Add T1 loss terms
            paramsOut.A = A;
        end
    end
        methods
        function [TRList,Mxy,Mz] = compile(self,M0,params)
            % COMPILE: runs the model based on some input parameters
            params = self.parseParams(params);
            A = params.A;
            volumeFractions = params.volumeFractions;
            FaList = params.FaList;
            TRList = params.TRList;
            [TRList, Mxy, Mz] = self.evaluate(TRList,FaList,M0,A,volumeFractions);
        end
        function [x,resultParams,allParams,resnorm,residual,exitflag,output,lambda,jacobian]...
                = fitData(self,params,guess,xdata,ydata,varargin)
            %FITDATA: Fits some set of guess parameters to input data
            %following the given model
            p = inputParser();
            p.addOptional('lb',[])
            p.addOptional('ub',[])
            p.parse(varargin{:})
            xNames = fieldnames(guess);
            j = 1;
            xIndex = cell(size(xNames));
            for i = 1:numel(xNames)
                iFits = ~isnan(guess.(xNames{i}));
                xIndex{i} = find(iFits==1);
                for k = 1:numel(xIndex{i})
                    x0(j) = guess.(xNames{i})(xIndex{i}(k));
                    j = j+1;
                end
            end
            tmpFlipAnlge = params.FaList(:,1);
            params = self.parseParams(params);
            Y0 = ydata(:,1)./sin(tmpFlipAnlge);
            fun = @(x,xdata)self.fitFunction(...
                params,x,xNames,xIndex,Y0);
            opts = params.fitOptions;
            [x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                lsqcurvefit(fun,x0,xdata,ydata,...
                [p.Results.lb],[p.Results.ub],opts);
            resultParams = guess;
            allParams = params;
            j = 1;
            for i = 1:numel(xNames)
                if strcmp(xNames{i},'FaList')
                    resultParams.(xNames{i}) = x(j);
                    allParams.(xNames{i}) = x(j);
                    j = j+1;
                else
                    for k = 1:numel(xIndex{i})
                        resultParams.(xNames{i})(xIndex{i}(k)) = x(j);
                        allParams.(xNames{i})(xIndex{i}(k)) = x(j);
                        j = j+1;
                    end
                end
            end
        end
        function DataCompare(self,params,xdata,ydata)
            % DATACOMPARE: a function for comparing some data with a set of
            % parameters and the model, the inital condition for the model
            % is taken from the data.
            M0 = self.splitM0(ydata(:,1),params);
            M0 = M0./sin(params.FaList(:,1));
            [TRList,Mxy,~] = self.compile(M0,params);
            legendVals = cell(size(Mxy,1),1);
            figure
            for i = 1:size(Mxy,1)
                tmpLine = plot(TRList,Mxy(i,:));
                hold on
                plot(xdata,ydata(i,:),'o','MarkerEdgeColor',tmpLine.Color);
                legendVals{2*i-1} = ['Fit Pool ',char(i+'A'-1)];
                legendVals{2*i} = ['Data Pool ',char(i+'A'-1)];
            end
            hold off
            xlabel('Time (sec)')
            ylabel('Signal (arb)')
            legend(legendVals)
        end
        function Y = fitFunction(self,params,x,xNames,xIndex,Y0)
            % fitFunction packs the parameter in params and x up and evaluates
            % using the evaluate funnction over some time (tSpan) with some
            % initial value (Y0)
            j = 1;
            for i = 1:numel(xNames)
                for k = 1:numel(xIndex{i})
                    params.(xNames{i})(xIndex{i}(k)) = x(j);
                    % Check if fitting flip angle (there mus be a better
                    % way to do this
                    if strcmp(xNames{i}, 'FaList')
                        params.(xNames{i}) =...
                            repmat(x(j),size(params.(xNames{i})));
                    end
                    j = j+1;
                end
            end
            % Split out the multiple Physical compartments
            Y0 = self.splitM0(Y0,params);
            [~, Y, ~] = self.compile(Y0,params);
        end
    end
    
    methods (Access = private)
        function [TRList, Mxy, Mz] = evaluate(~,TRList,FaList,M0,A,volumeFractions)
            % EVALUATE: runs the model based on some input parameters
            fun = @(t,y)A*y;
            tmpMz = zeros(size(FaList));
            tmpMxy = zeros(size(FaList));
            tmpMz(:,1) = M0.*cos(FaList(:,1));
            tmpMxy(:,1) = (M0.*sin(FaList(:,1)));
            for i = 2:length(TRList)
                % the transpose on Mz dose not matter as matlab
                % automatically converts the vector to match matrix
                % multiplication conventions but is done for clairity to
                % match with fun as it is declared above
                [~,Y] = ode45(fun,[TRList(i-1),TRList(i)],tmpMz(:,i-1));
                tmpMz(:,i) = Y(end,:).';
                tmpMxy(:,i) = sin(FaList(:,i)).*tmpMz(:,i);
                tmpMz(:,i) = cos(FaList(:,i)).*tmpMz(:,i);
            end
            nP = length(volumeFractions);
            nC = size(A,1)/nP;
            volumeFractions = reshape(repmat(volumeFractions,nC,1),nC*nP,1);
            tmpMxy = volumeFractions.*tmpMxy;
            tmpMz = volumeFractions.*tmpMz;
            Mxy = zeros(nC,length(TRList));
            Mz = zeros(nC,length(TRList));
            for i = 1:nC
                Mxy(i,:) = sum(tmpMxy(i:nC:end,:),1);
                Mz(i,:) = sum(tmpMz(i:nC:end,:),1);
            end
        end
        function M0out = splitM0(~,M0in,params)
            % Split out the multiple Physical compartmets. Right now it
            % will fill just the first compartment. A method to properly
            % fill all of the compartmetns with their respective fractions
            % of the input signal is unclear.
            nC = size(params.ExchangeTerms,2);
            nP = size(params.PerfusionTerms,2);
            M0out = [M0in;zeros((nC*(nP-1)),1)]./params.volumeFractions(1);
        end
    end
end