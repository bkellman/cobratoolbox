function [controlFluxs, objFlux] = nRobustnessAnalysis(model, controlRxns, nPoints, plotResFlag, objRxn,objType,emperical)
% Performs robustness analysis for a pair of reactions of
% interest and an objective of interest
%
% USAGE:
%
%    [controlFlux1, controlFlux2, objFlux] = doubleRobustnessAnalysis(model, controlRxn1, controlRxn2, emperical, nPoints, plotResFlag, objRxn, objType)
%
% INPUTS:
%    model:          COBRA model structure
%    controlRxns:    Vector of reactions of interest whose value is to be controlled
%
% OPTIONAL INPUTS:
%    emperical:      numerical vector of emperical growth, if specified, z-axis is -log(error) 
%    nPoints:        Number of flux values per dimension (Default = 20)
%    plotResFlag:    Plot results (Default = true)
%    objRxn:         Objective reaction to be maximized (Default = whatever
%                    is defined in model)
%    objType:        Maximize ('max') or minimize ('min') objective
%                    (Default = 'max')
%
% OUTPUTS:
%    controlFluxs:   Vector of flux values within the range of the maximum and minimum for
%                    reaction of interest
%    objFlux:        Optimal values of objective reaction at each control
%                    reaction flux value
%
% .. Authors: - Ben Kellman & Matt Schinn 6/13/19

if (nargin < 3)
    nPoints = 20;
end
if (nargin < 4)
    plotResFlag = true;
end
if (nargin > 5)
    baseModel = changeObjective(model,objRxn);
else
    baseModel = model;
end
if (nargin <6)
    objType = 'max';
end
if (nargin <7)
    emperical=[];
end

solMins = {};
solMaxs = {};

for rxnI=controlRxns
    if (findRxnIDs(model,rxnI))
        tmpModel = changeObjective(model,rxnI);
        solMins{rxnI} = optimizeCbModel(tmpModel,'min');
        solMaxs{rxnI} = optimizeCbModel(tmpModel,'max');
    else
        error('Control reaction 1 does not exist!');
    end
end

objFlux = [];
controlFluxs = [];
for rxnI=controlRxns
   controlFluxs = linspace(solMin{rxnI}.f,solMax{rxnI}.f,nPoints)'; 
end

showprogress(0,'Double robustness analysis in progress ...');
for i=1:nPoints
    for j = 1:nPoints
        showprogress(((i-1)*nPoints+j)/nPoints^2);
        modelControlled = changeRxnBounds(baseModel,controlRxn1,controlFlux1(i),'b');
        for rxnI=controlRxns[-1]
            modelControlled = changeRxnBounds(modelControlled,rxnI,controlFlux(rxnI)(j),'b');
        end
        solControlled = optimizeCbModel(modelControlled,objType);
        if length(emperical)>0
            resp = -log(solControlled.f - emperical)
        else
            resp = solControlled.f
        end
        objFlux(i,j) = resp;
    end
end


