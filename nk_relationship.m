% input nk_wav(波长序列，单位nm)，nk_int(强度序列，单位a.u.)，nk_flag(n,k,代表输入的是n还是k)
% 添加d_ome大小和，ome_num，分别表示ome序列的间隔和长度
function nk_out = nk_relationship(nk_wav,nk_int,nk_flag,d_ome,ome_num)    
    %% 定义频率，导入k值或者n值
    nk_ome = 1240./nk_wav; % unit nm
    %% 创建间隔为d_ome的序列，并进行插值，得到对应的强度序列
    ome_std = linspace(d_ome,ome_num*d_ome,ome_num);
    int_std = interp1(nk_ome,nk_int,ome_std,"nearest","extrap");
    int_std(int_std<0) = 0;
%     figure
%     plot(ome_std,int_std);
    %% 判断nk_flag是n还是k
    ome_0 = ome_std;
    for jj = 1:ome_num
        P = 0;
        for ii = 1:ome_num
            if ii == jj
            else
                if nk_flag == 'n'
                    % Σn(w).*dw./(w0.^2-w.^2)
                    cot = (int_std(ii).*d_ome./(ome_0(jj).^2-ome_std(ii).^2));
                elseif nk_flag == 'k'
                    % Σk(w).*w.*dw./(w.^2-w0.^2)
                    cot = (int_std(ii).*ome_std(ii).*d_ome./(ome_std(ii).^2-ome_0(jj).^2));
                else error('Wrong nk_flag');
                end
                P = P+cot;
            end
        end
        if nk_flag == 'n'
            out_temp(jj) = 2.*ome_0(jj)./pi.*P;
        else 
            out_temp(jj) = 1+2./pi.*P;
        end
    end
    nk_out = out_temp;
    figure 
    plot(ome_0,nk_out,'r-');
    title('transfer out');
end