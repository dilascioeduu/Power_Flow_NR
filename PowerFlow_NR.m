clc
clear all
close all 

% 1 - Entrada de dados das barras

% fprintf('=CONFIGURAÇÃO DAS BARRAS=\n');
% nb = input('Digite o número total de barras: '); 
% 
% barras = zeros(nb, 6);
% % P/ barras de carga (PQ), as potências especificadas precisam ser
% % negativas
% for i = 1:nb
%     fprintf('\n------ Barra %d ------\n', i);
%     barras(i,1) = i;
%     barras(i,2) = input('  Tipo (1-Slack, 2-PV, 3-PQ): ');
%     barras(i,3) = input('  Tensão (V) inicial [pu]: ');
%     barras(i,4) = input('  Ângulo (theta) inicial [graus]: ') * (pi/180);
%     barras(i,5) = input('  Potência Ativa (P) especificada [pu]: ');
%     barras(i,6) = input('  Potência Reativa (Q) especificada [pu]: ');
% end
% 
% % 2 - Entrada de dados das linhas
% 
% fprintf('\n=CONFIGURAÇÃO DAS LINHAS=\n');
% nl = input('Digite o número de linhas de transmissão: ');
% 
% linhas = zeros(nl, 4);
% for k = 1:nl
%     fprintf('\n------Linha %d------\n', k);
%     linhas(k,1) = input('  Barra de origem (De): ');
%     linhas(k,2) = input('  Barra de destino (Para): ');
%     linhas(k,3) = input('  Resistência R [pu]: ');
%     linhas(k,4) = input('  Reatância X [pu]: ');
% end


arquivo_excel = 'sistema.xlsx';

fprintf('=== IMPORTANDO DADOS DO EXCEL ===\n');

% Configuração das bases globais
Sbase = 100; % MVA
fprintf('=== PROCESSANDO DADOS BRUTOS NO PADRÃO IEEE ===\n');

% 1. Leitura direta das abas do Excel
tabela_raw_barras = readtable(arquivo_excel, 'Sheet', 'Barras');
tabela_raw_linhas = readtable(arquivo_excel, 'Sheet', 'Linhas');

nb = size(tabela_raw_barras, 1);
nl = size(tabela_raw_linhas, 1);

% 2. Pré-alocação da matriz final do algoritmo: [ID, Tipo, V, Theta, P_pu, Q_pu]
barras = zeros(nb, 6);

for i = 1:nb
    barras(i,1) = tabela_raw_barras.Bus_No(i); % Copia o ID da Barra
    
    % Tradução automática do código de operação IEEE para o Tipo que
    % estamos utilizando (1 - Slack; 2 - PV; 3 - PQ)
    ieee_code = tabela_raw_barras.Bus_Code(i);
    if ieee_code == 1
        barras(i,2) = 1; % 1 vira Slack 
    elseif ieee_code == 2
        barras(i,2) = 2; % 2 vira PV 
    else
        barras(i,2) = 3; % 0 vira PQ 
    end
    
    barras(i,3) = tabela_raw_barras.V_pu(i); % Módulo de tensão inicial
    barras(i,4) = 0; % Ângulo inicial (0 rad)
    
    % Adaptando as potências líquidas da tabela: (Geração - Carga) / Sbase
    barras(i,5) = (tabela_raw_barras.Gen_MW(i) - tabela_raw_barras.Load_MW(i)) / Sbase;
    barras(i,6) = (tabela_raw_barras.Gen_Mvar(i) - tabela_raw_barras.Load_Mvar(i)) / Sbase;
end

% 3. Montagem direta da matriz de linhas 
linhas = [tabela_raw_linhas.From_Bus, ...
          tabela_raw_linhas.To_Bus, ...
          tabela_raw_linhas.R_pu, ...
          tabela_raw_linhas.X_pu, ...
          tabela_raw_linhas.B_pu];

fprintf('>> Dados convertidos para pu e condicionados com sucesso!\n');

% Configurações de execução

tol = input('\nDigite a tolerância (ex: 1e-4): ');
max_iter = 20;
tic
% 3 - Formação da matriz de admitância (Y BUS)

Ybus = zeros(nb, nb);
for k = 1:nl
    i = find(barras(:,1) == linhas(k,1)); % Resgata as barras de origem (De) que o cabo k interliga
    j = find(barras(:,1) == linhas(k,2)); % Resgata as barras de destino (Para) que o cabo k interliga
    z = linhas(k,3) + 1j*linhas(k,4); % Monta a impedância (R + jX)
    B_shunt = linhas(k,5);
    y = 1/z; % Transforma em admitância (Y = 1/Z)
    Ybus(i,j) = -y; % Elementos fora da diagonal principal são negativos
    Ybus(j,i) = -y;
    Ybus(i,i) = Ybus(i,i) + y + 1j*B_shunt;
    Ybus(j,j) = Ybus(j,j) + y + 1j*B_shunt;
end
Ym = abs(Ybus); % Armazena o valor do módulo de cada elemento da matriz
Yth = angle(Ybus); % Armazena o valor do ângulo de cada elemento da matriz

% 4 - Início do processo iterativo (NR)

V = barras(:,3); % Cria vetores coluna das magnitudes das tensões e ângulos
theta = barras(:,4);
iter = 0;
conv = false;

while ~conv && iter < max_iter
    iter = iter + 1;
    Pcalc = zeros(nb, 1);
    Qcalc = zeros(nb, 1);

    for i = 1:nb
        for j = 1:nb
            ang = Yth(i,j) + theta(j) - theta(i); % alpha = gamaij + thetai - thetaj
            Pcalc(i) = Pcalc(i) + Ym(i,j) * V(i) * V(j) * cos(ang); % Potência ativa injetada na barra
            Qcalc(i) = Qcalc(i) - Ym(i,j) * V(i) * V(j) * sin(ang); % Potência reativa injetada na barra
        end
    end

    dP = barras(:,5) - Pcalc; % Vetores de Mismatches (resíduos) (dP = Pesp - Pcalc)
    dQ = barras(:,6) - Qcalc;
    idx_P = find(barras(:,2) ~= 1); % Procura na matriz onde a barra não é referência (PQ ou PV precisam de dP)
    idx_Q = find(barras(:,2) == 3); % Procura na matriz onde a barra É PQ (só ela precisa de dQ)
    mis = [dP(idx_P); dQ(idx_Q)]; % Cria o vetor de Mismatches deltaP e deltaQ

    if max(abs(mis)) < tol % Verifica o critério de parada
        conv = true;
        break;
    end

    % Montagem da Matriz Jacobiana
    H = zeros(nb, nb);
    N = zeros(nb, nb);
    M = zeros(nb, nb); 
    L = zeros(nb, nb);
    for i = 1:nb
        for j = 1:nb
            ang = Yth(i,j) + theta(j) - theta(i);
            if i ~= j % Verifica se estamos fora da diagonal principal
                H(i,j) = -V(i) * V(j) * Ym(i,j) * sin(ang);
                N(i,j) =  V(i) * Ym(i,j) * cos(ang);
                M(i,j) = -V(i) * V(j) * Ym(i,j) * cos(ang);
                L(i,j) = -V(i) * Ym(i,j) * sin(ang);
            else % Cálcuços utilizando as simplificações
                H(i,i) = -Qcalc(i) - (V(i)^2 * imag(Ybus(i,i)));
                N(i,i) = (Pcalc(i) / V(i)) + (V(i) * real(Ybus(i,i)));
                M(i,i) = Pcalc(i) - (V(i)^2 * real(Ybus(i,i)));
                L(i,i) = (Qcalc(i) / V(i)) - (V(i) * imag(Ybus(i,i)));
            end
        end
    end

    % Montagem da Jacobiana eliminando as equações da barra Slack e as
    % tensões das barras PV (fixas)
    J = [H(idx_P, idx_P), N(idx_P, idx_Q); M(idx_Q, idx_P), L(idx_Q, idx_Q)];
    
    [L_mat, U_mat] = lu(J); % Fatoração LU
    dx = U_mat \ (L_mat \ mis); % Faz a substituição direta (L*y = dMismatches) e 
    % a retroativa (U*dx=y) em uma só operação

    n_p = length(idx_P); 
    theta(idx_P) = theta(idx_P) + dx(1:n_p); % Atualização do estado angular para a próxima iteração
    V(idx_Q) = V(idx_Q) + dx(n_p+1:end); % Atualização do estado de tensão para a próxima iteração
end

% 5 - Cálculo dos fluxos (P e Q)
fprintf('\n=PROCESSAMENTO DE FLUXOS E PERDAS=\n');

P_fluxo = zeros(nl, 2); % [P_de_para, P_para_de]
Q_fluxo = zeros(nl, 2); % [Q_de_para, Q_para_de]
P_perda = zeros(nl, 1);
Q_perda = zeros(nl, 1);

for k = 1:nl
    i = find(barras(:,1) == linhas(k,1)); %i = linhas(k,1); 
    % Resgata as barras de origem (De) que o cabo k interliga
    j = find(barras(:,1) == linhas(k,2)); % j = linhas(k,2); 
    % Resgata as barras de destino (Para) que o cabo k interliga
    y_linha = 1 / (linhas(k,3) + 1j*linhas(k,4)); % Adamitância do próprio cabo
    b_metade = 1j * (linhas(k,5) / 2); % Metade da susceptância Shunt (modelo pi)
    
    Vi = V(i) * exp(1j*theta(i));
    Vj = V(j) * exp(1j*theta(j));
    
    % Potência complexa S = V * conj(I)
    S_ij = Vi * conj((Vi - Vj) * y_linha + Vi * b_metade);
    S_ji = Vj * conj((Vj - Vi) * y_linha + Vj * b_metade);
    
    % Separação em Ativa e Reativa
    P_fluxo(k,1) = real(S_ij);
    Q_fluxo(k,1) = imag(S_ij);
    P_fluxo(k,2) = real(S_ji);
    Q_fluxo(k,2) = imag(S_ji);
    
    % Perdas na linha
    P_perda(k) = real(S_ij + S_ji);
    Q_perda(k) = imag(S_ij + S_ji);
    
    Perdas_Totais_Ativas_pu = sum(P_perda);
    Perdas_Totais_Reativas_pu = sum(Q_perda);

    % Conversão para grandezas reais usando a Potência Base (MW e Mvar)
    Perdas_Totais_Ativas_MW = Perdas_Totais_Ativas_pu * Sbase;
    Perdas_Totais_Reativas_Mvar = Perdas_Totais_Reativas_pu * Sbase;
end

% Definição das bases para conversão de grandezas reais
% Sbase = input('\nDigite a potência base do sistema(MVA):'); % Potência Base do Sistema (MVA)
Vbase = input('\nDigite a tensão base do sistema(kV):'); % Tensão Base das Barras (kV)

% Pré-alocação das colunas de fluxo real
Nom_kV = ones(nb, 1) * Vbase;
PU_Volt = V;
Volt_kV = V * Vbase;
Angle_Deg = rad2deg(theta);
Gen_MW = zeros(nb, 1);
Gen_Mvar = zeros(nb, 1);
Load_MW = zeros(nb, 1);
Load_MVar = zeros(nb, 1);

From_Number = linhas(:,1);
To_Number = linhas(:,2);

MW_From = P_fluxo(:,1) * Sbase;
MVar_From = Q_fluxo(:,1) * Sbase;
MVA_From = sqrt(MW_From.^2 + MVar_From.^2);

MW_To = P_fluxo(:,2) * Sbase;
MVar_To = Q_fluxo(:,2) * Sbase;
MVA_To = sqrt(MW_To.^2 + MVar_To.^2);

MW_Loss = P_perda * Sbase;
MVar_Loss = Q_perda * Sbase;

Tabela_Linhas_PW = table(From_Number, To_Number, ...
    MW_From, MVar_From, MVA_From, MW_To, MVar_To, MVA_To, MW_Loss, MVar_Loss, ...
    'VariableNames', {'From', 'To', ...
                      'MW_ik', 'Mvar_ik', 'MVA_ik', 'MW_ki', 'MVar_ki', 'MVA_ki', 'MW_Loss', 'Mvar_Loss'});
                  

% Separação física dos fluxos de acordo com o tipo de barra
for i = 1:nb
    tipo = barras(i,2);
    
    if tipo == 1 % Barra Slack (Referência)
        Gen_MW(i) = Pcalc(i) * Sbase;
        Gen_Mvar(i) = Qcalc(i) * Sbase;
        
    elseif tipo == 2 % Barra PV (Gerador)
        Gen_MW(i) = barras(i,5) * Sbase; % P especificado original
        Gen_Mvar(i) = Qcalc(i) * Sbase;     % Q calculado pela rede
        
    elseif tipo == 3 % Barra PQ (Carga)
        % Inverte o sinal (-) do MATLAB para virar consumo (+) no PowerWorld
        Load_MW(i) = -barras(i,5) * Sbase; 
        Load_MVar(i) = -barras(i,6) * Sbase;
    end
end

% 6. APRESENTAÇÃO DOS RESULTADOS
fprintf('\nIterações para convergência: %d\n', iter);
toc

% 1. Puxa os dados estáticos originais diretamente do Excel bruto
Gen_MW = tabela_raw_barras.Gen_MW;
Gen_Mvar = tabela_raw_barras.Gen_Mvar;
Load_MW = tabela_raw_barras.Load_MW;
Load_Mvar = tabela_raw_barras.Load_Mvar;

% 2. Reconstitui APENAS as variáveis calculadas pelo Newton-Raphson
 type = cell(nb,1);

for i = 1:nb
    tipo = barras(i,2); % Identifica o tipo de barra (1=Slack, 2=PV, 3=PQ)
    
    if tipo == 1 
        % Barra Slack: Calcula a geração real ativa e reativa necessárias
        Gen_MW(i)   = (Pcalc(i) * Sbase) + Load_MW(i);
        Gen_Mvar(i) = (Qcalc(i) * Sbase) + Load_Mvar(i);
        type{i} = 'SLACK';
        

    elseif tipo == 2 
        % Barra PV: Mantém o Gen_MW fixo do Excel e calcula o Gen_Mvar gerado
        Gen_Mvar(i) = (Qcalc(i) * Sbase) + Load_Mvar(i);
        type{i} = 'PV';
        
    elseif tipo == 3
        % Barra PQ: Não gera nada, garante que as colunas de geração fiquem zeradas
        Gen_MW(i) = 0;
        Gen_Mvar(i) = 0;
        type{i} = 'PQ';
       
    end
end

% 3. Preparação das colunas de tensão e ângulo
PU_Volt = V;
Volt_kV = V * 138; % Altere 138 para a tensão base do seu sistema se necessário
Angle_Deg = rad2deg(theta);  

% 4. Montagem da tabela final idêntica ao PowerWorld
Tabela_PowerWorld = table(tabela_raw_barras.Bus_No, type, PU_Volt, Volt_kV, Angle_Deg, ...
    Gen_MW, Gen_Mvar, Load_MW, Load_Mvar, ...
    'VariableNames', {'Name', 'Type', 'PU_Volt', 'Volt_kV', 'Angle_Deg', 'Gen_MW', 'Gen_Mvar', 'Load_MW', 'Load_Mvar'});

fprintf('\n=====================================================================================================================\n');
fprintf('                                                   BUSES \n');
fprintf('=====================================================================================================================\n');
disp(Tabela_PowerWorld);

% Tabela_Barras = table((1:nb)', V, rad2deg(theta), Pcalc, Qcalc, ...
%     'VariableNames', {'Barra', 'V_pu', 'Angulo_deg', 'P_Gerada_pu', 'Q_Gerada_pu'});
% disp('--- ESTADO DAS BARRAS ---');
% disp(Tabela_Barras);
% 
% Tabela_Fluxos = table(linhas(:,1), linhas(:,2), P_fluxo(:,1), Q_fluxo(:,1), P_perda, Q_perda, ...
%     'VariableNames', {'De', 'Para', 'P_Fluxo_pu', 'Q_Fluxo_pu', 'Perda_P_pu', 'Perda_Q_pu'});
% disp('--- FLUXOS E PERDAS NAS LINHAS ---');
% disp(Tabela_Fluxos);

% Tabela_Estilo_PW = table((1:nb)', Nom_kV, PU_Volt, Volt_kV, Angle_Deg, ...
%     Gen_MW, Gen_Mvar, Load_MW, Load_MVar, ...
%     'VariableNames', {'Name', 'Nom_kV', 'PU_Volt', 'Volt_kV', 'Angle_Deg', ...
%                       'Gen_MW', 'Gen_Mvar', 'Load_MW', 'Load_Mvar'});
% 
% fprintf('\n===================================================================\n');
% fprintf('         RESULTADOS FORMATADOS EM GRANDEZA REAL (ESTILO POWERWORLD) \n');
% fprintf('===================================================================\n');
% disp(Tabela_Estilo_PW);

fprintf('\n=====================================================================================================================\n');
fprintf('                                                 POWER FLOW \n');
fprintf('\n=====================================================================================================================\n');
disp(Tabela_Linhas_PW);

% Imprime o resumo global de perdas no terminal
fprintf('\n=====================================================================================================================\n');
fprintf('                         RESUMO GLOBAL DE PERDAS DA REDE                        \n');
fprintf('=====================================================================================================================\n');
fprintf('Perdas Ativas Totais (P_loss): %.4f pu (%.2f MW)\n', ...
    Perdas_Totais_Ativas_pu, Perdas_Totais_Ativas_MW);
fprintf('Perdas Reativas Totais (Q_loss): %.4f pu (%.2f Mvar)\n', ...
    Perdas_Totais_Reativas_pu, Perdas_Totais_Reativas_Mvar);
fprintf('=====================================================================================================================\n');