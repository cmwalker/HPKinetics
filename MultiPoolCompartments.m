classdef MultiPoolCompartments < HPKinetics.MultiPool
    %MULTIPOOLCOMPARTMENT Multi-Physical and Chemical Pool models
    %   This model builds off the MultiPool model for multiple chemical
    %   pools but expands it for multiple physical compartments.
    
    properties
    end
    
    methods
        function defaults(self)
            % DEFAULTS explains the default values for each parameter
            names = {'PerfusionTerms','volumeFractions'};
            discriptions = {' A Row Vector of perfusion Exchange Constnats for each chemical pool.'...
                ' A Row Vector of volme fraction for each physical pool. Only one value can be use if all pools have the same volume fraction.'};
            defaultsVals = {'0','1'};
            fprintf('*Note* all terms must be a vector of size 1 x N where N is the number of chemical Pools\n')
            for i = 1:numel(names)
                fprintf('''%s'': %s\n Default Vaule: %s\n',...
                    names{i},discriptions{i},defaultsVals{i});
            end
            defaults@HPKinetics.MultiPool(self)
        end
        function paramsOut = parseParams(self,paramsIn)
            % parseParams: a function to fill default param values if they are
            % not defined
            default = struct('PerfusionTerms',0,'volumeFractions',1);
            tmpNames = fieldnames(default);
            paramsOut = paramsIn;
            for i = 1:numel(tmpNames)
                if ~isfield(paramsOut,tmpNames{i})
                    paramsOut.(tmpNames{i}) = default.(tmpNames{i});
                end
            end
            paramsOut = parseParams@HPKinetics.MultiPool(self,paramsOut);
            N = size(paramsOut.ExchangeTerms,1);
            kve = paramsOut.PerfusionTerms;
            if(length(paramsOut.volumeFractions)==size(paramsOut.ExchangeTerms,1))
                % if the corect number of ve terms are passed in continue
                ve = zeros(N,1)+paramsOut.volumeFractions;
                if sum(ve)>1
                    warning('The total volume of the system is %f, which is greater than unity\n',sum(ve))
                end
            elseif (length(paramsOut.volumeFractions)==(size(paramsOut.ExchangeTerms,1)-1))
                % If not enough volume fractions are passed in fill the
                % last on with the remainder
                ve = paramsOut.volumeFractions;
                ve(end+1) = 1-sum(paramsOut.volumeFractions);
            else
                error('Not enough volume fractions passed in %d physical pools, %d volume fractions\n',...
                    length(paramsOut.volumeFractions),size(paramsOut.ExchangeTerms,1))
            end
            A = paramsOut.A;
            for i = 1:N
                for j = 1:N
                    if(i == j)
                        A(i,i) = A(i,i)-kve(i)/ve(i);
                    end
                end
            end
            paramsOut.A = A;
            paramsOut.kve = paramsOut.PerfusionTerms;
            paramsOut.ve = paramsOut.volumeFractions;
        end
        function [TRList,Mxy,Mz] = compile(self,M0,params)
            % EVALUATE: runs the model based on some input parameters
            params = self.parseParams(params);
            A = params.A;
            FaList = params.FaList;
            TRList = params.TRList;
            [TRList, Mxy, Mz] = self.evaluate(...
                TRList,FaList,M0,A,params);
        end
        function DataCompare(self,params,xdata,ydata)
            ydata(1,:) = ydata(1,:)./params.volumeFractions(1);
            DataCompare@HPKinetics.MultiPool(self,params,xdata,ydata);
        end
        function [x,resultParams,allParams,resnorm,residual,exitflag,output,lambda,jacobian]...
                = fitData(self,params,guess,xdata,ydata,varargin)
            ydata(1,:) = ydata(1,:)./params.volumeFractions(1);
            [x,resultParams,allParams,resnorm,residual,exitflag,output,lambda,jacobian]...
                = fitData@HPKinetics.MultiPool(self,params,guess,xdata,ydata,varargin{:});
        end
    end
     methods (Access = private)
         function [TRList, Mxy, Mz] = evaluate(self,TRList,FaList,M0,A,params)
            % EVALUATE: runs the model based on some input parameters
            % thanks to matlabs matrix ordering for ODE45 kve have to
            % be transposed in the fun decleration
            fun = @(t,y)A*y;
            Mz = zeros(size(FaList));
            Mxy = zeros(size(FaList));
            Mz(1,:) = M0.*cos(FaList(1,:));
            Mxy(1,:) = params.ve.*M0.*sin(FaList(1,:));
            for i = 2:length(TRList)
                [~,Y] = ode45(fun,[TRList(i-1),TRList(i)],Mz(i-1,:).');
                Mz(i,:) = Y(end,:);
                Mxy(i,:) = sin(FaList(i,:)).*(params.ve.*Mz(i,:));
                Mz(i,:) = cos(FaList(i,:)).*Mz(i,:);
            end
        end
     end
end

