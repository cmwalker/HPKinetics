classdef MultiPoolTofftsGammaVIF < HPKinetics.MultiPoolToffts
    %MULTIPOOLTOFFTSGAMMAVIF A chemical exchange model assuming two pooled
    %Tofts model of perfusion
    %   parameters Values
    %* ExchangeTerms - A Matrix defining chemical Exchange. Defalt: 0
    %* T1s - A row vector of T1 decay terms. Default: 100
    %* FaList - A matrix of excitation angles. Default: 0
    %* TRList - A matrix of excitation times. Default: 0
    %* t0 - A row vector for delivery delay of each metabolite. Default: 0
    %* gammaPdfA - A row vector for shape term alpha of each metabolite. Default: 2.8
    %* gammaPdfB - A row vector for shape term beta of each metabolite. Default: 4.5
    %* ScaleFactor - A row vector for each metabolite's VIF scale factor. Default: 1
    %* fitOptions - A matlab fit option structure. Default: optimset(''lsqcurvefit'')
    %* PerfusionTerms - A row vector for each metabolite's extravisation rate. Default: 0
    %* volumeFractions - A row vector for each metabolite's volume fraction. Default: 1
    %   There is NO imput validation for the parameters passed in, for more
    %   detail on the assumed data structur of these parameters use the
    %   defaults function
    properties
    end
    methods
        function defaults(self)
            % DEFAULTS explains the default values for each parameter
            names = {'t0','gammaPdfA','gammaPdfB','scaleFactor'};
            discriptions = {'A  Row vector of time delays for each metabolite'...
                ' A  Row vector of shape term Alpha, set this to zero to have no VIF for a chemical pool'...
                ' A  Row vector of shape term Beta, this cannot be zero and will be set to 1e-40 if zero is used'...
                ' A  Row vector of Scale Factor to be applied to the VIF'};
            defaultsVals = {'0','2.8','4.5','1'};
            fprintf('*Note* all terms must be a vector of size 1 x N where N is the number of chemical Pools\n')
            for i = 1:numel(names)
                fprintf('''%s'': %s\n Default Vaule: %s\n',...
                    names{i},discriptions{i},defaultsVals{i});
            end
            defaults@HPKinetics.MultiPoolToffts(self);
        end
        function paramsOut = parseParams(self,paramsIn)
            % parseParams: Parses the input shape terms of a gamma variate
            % for the VIF, each term should be a vector with shape terms
            % for each chemical species
            
            % Fill Default Values
            default = struct('t0',0,'gammaPdfA',2.8,...
                'gammaPdfB',4.5,'scaleFactor',1);
            tmpNames = fieldnames(default);
            paramsOut = paramsIn;
            for i = 1:numel(tmpNames)
                if ~isfield(paramsOut,tmpNames{i})
                    paramsOut.(tmpNames{i}) = default.(tmpNames{i});
                end
            end
            % Build VIF
            paramsOut.VIF = @(t)paramsIn.scaleFactor.*...
                gampdf(t-paramsIn.t0,paramsIn.gammaPdfA,paramsIn.gammaPdfB);
            % Fill in parent Class defaults
            paramsOut = parseParams@HPKinetics.MultiPoolToffts(self,paramsOut);
        end
        function [TRList,Mxy,Mz] = compile(self,M0,params)
        % EVALUATE: runs the model based on some input parameters
        %params = self.parseParams(params);
        [TRList,Mxy,Mz] = compile@HPKinetics.MultiPoolToffts(self,M0,params);
        end
%         function [x,resultParams,allParams,resnorm,residual,exitflag,output,lambda,jacobian]...
%                 = fitData(self,params,guess,xdata,ydata,varargin)
%             p = inputParser();
%             p.addOptional('lb',[])
%             p.addOptional('ub',[])
%             p.parse(varargin{:})
%             xNames = fieldnames(guess);
%             j = 1;
%             xIndex = cell(size(xNames));
%             for i = 1:numel(xNames)
%                 iFits = ~isnan(guess.(xNames{i}));
%                 xIndex{i} = find(iFits==1);
%                 for k = 1:numel(xIndex{i})
%                     x0(j) = guess.(xNames{i})(xIndex{i}(k));
%                     j = j+1;
%                 end
%             end
%             params = self.parseParams(params);
%             Y0 = ydata(:,1)./sin(params.FaList(:,1));
%             fun = @(x,xdata)self.fitFunction(...
%                 params,x,xNames,xIndex,xdata,Y0);
%             opts = params.fitOptions;
%             [x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
%                 lsqcurvefit(fun,x0,xdata,ydata,...
%                 [p.Results.lb],[p.Results.ub],opts);
%             resultParams = guess;
%             allParams = params;
%             j = 1;
%             for i = 1:numel(xNames)
%                 for k = 1:numel(xIndex{i})
%                     resultParams.(xNames{i})(xIndex{i}(k)) = x(j);
%                     allParams.(xNames{i})(xIndex{i}(k)) = x(j);
%                     j = j+1;
%                 end
%             end
%         end
    end
    methods (Access = private)
%         function Y = fitFunction(params,x,xNames,xIndex,tSpan,Y0)
%             % fitFunction packs the parameter in params and x up and evaluates
%             % using the evaluate funnction over some time (tSpan) with some
%             % initial value (Y0)
%             j = 1;
%             for i = 1:numel(xNames)
%                 for k = 1:numel(xIndex{i})
%                     params.(xNames{i})(xIndex{i}(k)) = x(j);
%                     j = j+1;
%                 end
%             end
%             params.TRList = tSpan;
%             [~, Y, ~] = HypWright.Models.MultiPoolTofftsGammaVIF.compile(Y0,params);
%         end
    end
    
end
    