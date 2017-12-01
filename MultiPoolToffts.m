classdef MultiPoolToffts < HPKinetics.MultiPool
    %MULTIPOOLTOFFTS Multi-Physical and Chemical Pool model following Tofts
    %   This model builds off the MultiPool model for multiple chemical
    %   pools but expands it to follow Tofts model for multiple physical
    %   compartments.
    
    properties
    end
    
    methods
        function defaults(self)
            % DEFAULTS explains the default values for each parameter
            names = {'PerfusionTerms','volumeFractions','VIF'};
            discriptions = {' A Row Vector of perfusion Exchange Constnats for each chemical pool.'...
                ' A Row Vector of volme fraction for each chemical pool. Only one value can be use if all pools have the same volume fraction.'...
                ' A function of a time variable (t) in seconds that returns a Row vector for the VIF of each chemical pool at the time t.'};
            defaultsVals = {'0','1','@(t)0'};
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
            default = struct('PerfusionTerms',0,'volumeFractions',1,...
                'VIF',@(t)0);
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
            if(length(paramsOut.volumeFractions)==1)
                ve = zeros(N,1)+paramsOut.volumeFractions;
            else
                ve = paramsOut.volumeFractions;
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
            paramsOut.b = paramsOut.VIF;
            paramsOut.kve = paramsOut.PerfusionTerms;
            paramsOut.ve = paramsOut.volumeFractions;
        end
        function [TRList,Mxy,Mz] = compile(self,M0,params)
            % EVALUATE: runs the model based on some input parameters
            params = self.parseParams(params);
            A = params.A;
            b = params.b;
            FaList = params.FaList;
            TRList = params.TRList;
            [TRList, Mxy, Mz] = self.evaluate(...
                TRList,FaList,M0,A,b,params);
        end
    end
     methods (Access = private)
         function [TRList, Mxy, Mz] = evaluate(self,TRList,FaList,M0,A,b,params)
            % EVALUATE: runs the model based on some input parameters
            kve = params.kve;
            ve = params.ve;
            % thanks to matlabs matrix ordering for ODE45 kve and b have to
            % be transposed in the fun decleration
            fun = @(t,y)A*y+(kve.'/ve).*b(t).';
            Mz = zeros(size(FaList));
            Mxy = zeros(size(FaList));
            Mz(1,:) = M0.*cos(FaList(1,:));
            Mxy(1,:) = (params.ve*M0+(1-params.ve)*params.b(TRList(1))).*sin(FaList(1,:));
            for i = 2:length(TRList)
                [~,Y] = ode45(fun,[TRList(i-1),TRList(i)],Mz(i-1,:).');
                Mz(i,:) = Y(end,:);
                Mxy(i,:) = sin(FaList(i,:)).*(params.ve.*Mz(i,:)+...
                    (1-params.ve).*b(TRList(i)));
                Mz(i,:) = cos(FaList(i,:)).*Mz(i,:);
            end
        end
     end
end

