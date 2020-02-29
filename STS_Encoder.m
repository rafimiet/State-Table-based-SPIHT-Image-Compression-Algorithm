function Stream = STS_Encoder(S,rate,L)
    bitcount = 0;I = 0;size_LB = 0;size_UB = 0;
    while bitcount < rate*size(S,1)*size(S,2)
        if I == 0 % INITIALIZATION
            % wavelet decomposition level can be defined by users manually.
            [R,C] = size(S);
            Sig = zeros(R,C);
            BL = zeros(R/2,C/2);
            [br,bc] = size(BL);
            BL(br/2+1:br,1:bc) = 3;
            BL(1:br/2,bc/2+1:bc) = 3;
            % BL can have values of 0,3,4,5,7
            M = max(max(S));
            n = floor(log(M)/log(2));
            T = 2^n; % THRESHOLD
            Stream = {};
            iter = 1;
            BL(1:R/(2^(L)),1:C/(2^(L))) = 4;bitcount = 0;
            I = 1; mm = 1;
        elseif I == 1 % REFINEMENT PASS
            if bitcount >= rate*size(S,1)*size(S,2)
                break
            end
            [x,y] = find(BL >= 4); % Block by block and like this (1,1) (2,1)...(1,2),(2,2)...(1,3),(2,3)...
            n = numel(x); % Number of significant blocks
            BS = [];
            for ii= 1:n
                if bitcount >= rate*size(S,1)*size(S,2)
                    break
                end
                if Sig(2*x(ii)-1,2*y(ii)-1) == 0
                    if S(2*x(ii)-1,2*y(ii)-1) >= T
                        BS = [BS 1 0];
                        bitcount = bitcount + 2;
                        S(2*x(ii)-1,2*y(ii)-1) = S(2*x(ii)-1,2*y(ii)-1) - T;
                        Sig(2*x(ii)-1,2*y(ii)-1) = 1;
                    elseif S(2*x(ii)-1,2*y(ii)-1) <= -T
                        BS = [BS 1 1];
                        bitcount = bitcount + 2;
                        S(2*x(ii)-1,2*y(ii)-1) = S(2*x(ii)-1,2*y(ii)-1) + T;
                        Sig(2*x(ii)-1,2*y(ii)-1) = 1;
                    else
                        BS = [BS 0];
                        bitcount = bitcount + 1;
                    end
                else
                    if S(2*x(ii)-1,2*y(ii)-1) >= T
                        S(2*x(ii)-1,2*y(ii)-1) = S(2*x(ii)-1,2*y(ii)-1) - T;
                        BS = [BS 1];
                        bitcount = bitcount + 1;
                    elseif S(2*x(ii)-1,2*y(ii)-1) <= -T
                        BS = [BS 1];
                        bitcount = bitcount + 1;
                        S(2*x(ii)-1,2*y(ii)-1) = S(2*x(ii)-1,2*y(ii)-1) + T;
                    else
                        BS = [BS 0];
                        bitcount = bitcount + 1;
                    end
                end

                if Sig(2*x(ii)-1,2*y(ii)) == 0
                    if S(2*x(ii)-1,2*y(ii)) >= T
                        Sig(2*x(ii)-1,2*y(ii)) = 1;
                        S(2*x(ii)-1,2*y(ii)) = S(2*x(ii)-1,2*y(ii)) - T;
                        BS = [BS 1 0];
                        bitcount = bitcount + 2;
                    elseif S(2*x(ii)-1,2*y(ii)) <= -T
                        BS = [BS 1 1];
                        bitcount = bitcount + 2;
                        Sig(2*x(ii)-1,2*y(ii)) = 1;
                        S(2*x(ii)-1,2*y(ii)) = S(2*x(ii)-1,2*y(ii)) + T;
                    else
                        BS = [BS 0];
                        bitcount = bitcount + 1;
                    end
                else
                    if S(2*x(ii)-1,2*y(ii)) >= T
                        BS = [BS 1];
                        bitcount = bitcount + 1;
                        S(2*x(ii)-1,2*y(ii)) = S(2*x(ii)-1,2*y(ii)) - T;
                    elseif S(2*x(ii)-1,2*y(ii)) <= -T
                        BS = [BS 1];
                        bitcount = bitcount + 1;
                        S(2*x(ii)-1,2*y(ii)) = S(2*x(ii)-1,2*y(ii)) + T;
                    else
                        BS = [BS 0];
                        bitcount = bitcount + 1;
                    end
                end

                if Sig(2*x(ii),2*y(ii)-1) == 0
                    if S(2*x(ii),2*y(ii)-1) >= T
                        BS = [BS 1 0];
                        bitcount = bitcount + 2;
                        Sig(2*x(ii),2*y(ii)-1) = 1;
                        S(2*x(ii),2*y(ii)-1) = S(2*x(ii),2*y(ii)-1) - T;
                    elseif S(2*x(ii),2*y(ii)-1) <= -T
                        BS = [BS 1 1];
                        bitcount = bitcount + 2;
                        S(2*x(ii),2*y(ii)-1) = S(2*x(ii),2*y(ii)-1) + T;
                        Sig(2*x(ii),2*y(ii)-1) = 1;
                    else
                        BS = [BS 0];
                        bitcount = bitcount + 1;
                    end
                else
                    if S(2*x(ii),2*y(ii)-1) >= T
                        BS = [BS 1];
                        bitcount = bitcount + 1;
                        S(2*x(ii),2*y(ii)-1) = S(2*x(ii),2*y(ii)-1) - T;
                    elseif S(2*x(ii),2*y(ii)-1) <= -T
                        BS = [BS 1];
                        bitcount = bitcount + 1;
                        S(2*x(ii),2*y(ii)-1) = S(2*x(ii),2*y(ii)-1) + T;
                    else
                        BS = [BS 0];
                        bitcount = bitcount + 1;
                    end
                end

                if Sig(2*x(ii),2*y(ii)) == 0
                    if S(2*x(ii),2*y(ii)) >= T
                        BS = [BS 1 0];
                        bitcount = bitcount + 2;
                        Sig(2*x(ii),2*y(ii)) = 1;
                        S(2*x(ii),2*y(ii)) = S(2*x(ii),2*y(ii)) - T;
                    elseif S(2*x(ii),2*y(ii)) <= -T
                        BS = [BS 1 1];
                        bitcount = bitcount + 2;
                        Sig(2*x(ii),2*y(ii)) = 1;
                        S(2*x(ii),2*y(ii)) = S(2*x(ii),2*y(ii)) + T;
                    else
                        BS = [BS 0];
                        bitcount = bitcount + 1;
                    end
                else
                    if S(2*x(ii),2*y(ii)) >= T
                        BS = [BS 1];
                        bitcount = bitcount + 1;
                        S(2*x(ii),2*y(ii)) = S(2*x(ii),2*y(ii)) - T;
                    elseif S(2*x(ii),2*y(ii)) <= -T
                        BS = [BS 1];
                        bitcount = bitcount + 1;
                        S(2*x(ii),2*y(ii)) = S(2*x(ii),2*y(ii)) + T;
                    else
                        BS = [BS 0];
                        bitcount = bitcount + 1;
                    end
                end;
            end
            I = 2;k = 0;%whos BS;
            i = 1; j = 1;sub = 1;UB = [];
        elseif I == 2 % SIGNIFICANCE PASS
            while (i <= R/(2^(L+1)) && j <= C/(2^(L+1)))
                if bitcount >= rate*size(S,1)*size(S,2)
                    break
                end
                %%%%%%%%%%% INITIALIZE LB and UB %%%%%%%%%%%%%%%%%%%%
                if sub == 1 %LH
                    BL_addr_x = i; BL_addr_y = j + C/(2^(L+1)); % LH
                elseif sub == 2 %HL
                    BL_addr_x = i + R/(2^(L+1)); BL_addr_y = j; % HL
                else %HH
                    BL_addr_x = i + R/(2^(L+1)); BL_addr_y = j + C/(2^(L+1)); % HH
                end
                % BL_addr gives the address of blocks in LH5, HL5 and HH5
                % which are themselves already significant, as initialized
                if BL(BL_addr_x,BL_addr_y) == 7
                    % Can be either 7 or 4 or 5 coz initially set to 4
                    if sub < 3
                        sub = sub + 1;
                    else
                        sub = 1;
                        if j < C/(2^(L+1))
                            j = j + 1;
                        else
                            j = 1;
                            i = i + 1;
                        end
                    end
                else
                    if BL(BL_addr_x,BL_addr_y) == 4
                        [flag,LB] = descendant_map(BL_addr_x,BL_addr_y,S,T);
                        if flag == 0 % i.e it has insignificant descendants
                            BS = [BS 0 0 0 0];LB = []; % WHY four zeros if all blocks insig?
                            bitcount = bitcount + 4;
                            if sub < 3
                                sub = sub + 1;
                            else
                                sub = 1;
                                if j < C/(2^(L+1))
                                    j = j + 1;
                                else
                                    j = 1;
                                    i = i + 1;
                                end
                            end
                        else
                            UB = [BL_addr_x, BL_addr_y];
                            BL(BL_addr_x,BL_addr_y) = 5;
                        end
                    elseif BL(BL_addr_x,BL_addr_y) == 5
                        UB = [BL_addr_x, BL_addr_y];
                        LB = [];
                        for i3 = 1:4
                            if i3 == 1
                                x1 = 2*BL_addr_x-1;y1 = 2*BL_addr_y-1;
                            elseif i3 == 2
                                x1 = 2*BL_addr_x-1;y1 = 2*BL_addr_y;
                            elseif i3 == 3
                                x1 = 2*BL_addr_x;y1 = 2*BL_addr_y-1;
                            else
                                x1 = 2*BL_addr_x;y1 = 2*BL_addr_y;
                            end
                            if BL(x1,y1) ~= 7
                                LB = [LB;x1, y1];
                            end
                        end
                    end
                end
                %%%%%%%% LB and UB INITIALIZED %%%%%%%%%%%%%%%%%
                %%%%%%%% Generation of Bitstream and partial updation of BL %%%%%%
                while size(LB,1) ~= 0
                    if bitcount >= rate*size(S,1)*size(S,2)
                        break
                    end
%                     if size_LB < size(LB,1)
%                         size_LB = size(LB,1)
%                     end
                    % Temporary Bitstreams Generation
                    if BL(LB(1,1),LB(1,2)) <= 3
                        %if current block insignificant
                        BS1 = [];
                        for k = 1:4
                            if k == 1
                                x1 = 2*LB(1,1) - 1; y1 = 2*LB(1,2) - 1;
                            elseif k == 2
                                x1 = 2*LB(1,1) - 1; y1 = 2*LB(1,2);
                            elseif k == 3
                                x1 = 2*LB(1,1); y1 = 2*LB(1,2) - 1;
                            elseif k == 4
                                x1 = 2*LB(1,1); y1 = 2*LB(1,2);
                            end
                            if S(x1,y1) >= T
                                BS1 = [BS1 1 0];
                                Sig(x1,y1) = 1;
                                S(x1,y1) = S(x1,y1) - T;
                            elseif S(x1,y1) <= -T
                                BS1 = [BS1 1 1];
                                Sig(x1,y1) = 1;
                                S(x1,y1) = S(x1,y1) + T;
                            else
                                BS1 = [BS1 0];
                            end
                        end
                        if BL(LB(1,1),LB(1,2)) == 0
                            [flag,lb1] = descendant_map(LB(1,1),LB(1,2),S,T);
                            if max(BS1) == 1
                                BS = [BS 1 BS1];
                                bitcount = bitcount + 1 + numel(BS1);
                                BL(LB(1,1),LB(1,2)) = BL(LB(1,1),LB(1,2)) + 4;
                                if flag == 1
                                    BS = [BS 1];
                                    bitcount = bitcount + 1;
                                    BL(LB(1,1),LB(1,2)) = BL(LB(1,1),LB(1,2)) + 1;
                                    UB = [UB; LB(1,1),LB(1,2)];
                                    if size(LB,1) == 1
                                        LB = [];
                                    else
                                        LB = LB(2:end,1:end);
                                    end;
                                    LB = [lb1; LB];
                                else
                                    BS = [BS 0];lb1 = [];
                                    bitcount = bitcount + 1;
                                    if size(LB,1) == 1
                                        LB = [];
                                    else
                                        LB = LB(2:end,1:end);
                                    end;
                                end
                            else
                                if flag == 0
                                    BS = [BS 0];lb1 = [];
                                    bitcount = bitcount + 1;
                                    if size(LB,1) == 1
                                        LB = [];
                                    else
                                        LB = LB(2:end,1:end);
                                    end;
                                else
                                    BS = [BS 1 0 0 0 0 1];
                                    bitcount = bitcount + 6;
                                    BL(LB(1,1),LB(1,2)) = 5;
                                    UB = [UB; LB(1,1),LB(1,2)];
                                    if size(LB,1) == 1
                                        LB = [];
                                    else
                                        LB = LB(2:end,1:end);
                                    end;
                                    LB = [lb1; LB];
                                end
                            end
                        else % means == 3
                            if max(BS1) == 1
                                BS = [BS 1 BS1];
                                bitcount = bitcount + 1 + numel(BS1);
                                BL(LB(1,1),LB(1,2)) = BL(LB(1,1),LB(1,2)) + 4;
                            else
                                BS = [BS 0];
                                bitcount = bitcount + 1;
                            end
                            if size(LB,1) == 1
                                LB = [];
                            else
                                LB = LB(2:end,1:end);
                            end;
                        end
                    elseif BL(LB(1,1),LB(1,2)) == 4
                        [flag,lb1] = descendant_map(LB(1,1),LB(1,2),S,T);
                        if flag == 1
                            BS = [BS 1];
                            bitcount = bitcount + 1;
                            UB = [UB; LB(1,1),LB(1,2)];
                            BL(LB(1,1),LB(1,2)) = 5;
                            if size(LB,1) == 1
                                LB = [];
                            else
                                LB = LB(2:end,1:end);
                            end;
                            LB = [lb1; LB];
                        else
                            BS = [BS 0];lb1 = [];
                            bitcount = bitcount + 1;
                            if size(LB,1) == 1
                                LB = [];
                            else
                                LB = LB(2:end,1:end);
                            end;
                        end
                    elseif BL(LB(1,1),LB(1,2)) == 5
                        LB1 = [];
                        for i3 = 1:4
                            if i3 == 1
                                x1 = 2*LB(1,1)-1;y1 = 2*LB(1,2)-1;
                            elseif i3 == 2
                                x1 = 2*LB(1,1)-1;y1 = 2*LB(1,2);
                            elseif i3 == 3
                                x1 = 2*LB(1,1);y1 = 2*LB(1,2)-1;
                            else
                                x1 = 2*LB(1,1);y1 = 2*LB(1,2);
                            end
                            if BL(x1,y1) ~= 7
                                LB1 = [LB1;x1, y1];
                            end
                        end
                        UB = [UB; LB(1,1),LB(1,2)];
                        if size(LB,1) == 1
                            LB = [];
                        else
                            LB = LB(2:end,1:end);
                        end;
                        LB = [LB1;LB];
                    else % i.e. = 7
                        if size(LB,1) == 1
                            LB = [];
                        else
                            LB = LB(2:end,1:end);
                        end;
                    end
                end
                %%%%%%% Final updation of BL using UB %%%%%%%%%%%
                while size(UB,1) ~= 0
                    n = size(UB,1);
%                     if size_UB < size(UB,1)
%                         size_UB = size(UB,1)
%                     end
                    xu = UB(n,1);   yu = UB(n,2);
                    des = 1;
                    for i1 = 1:4
                        if i1 == 1
                            xx = 2*xu-1;    yy = 2*yu-1;
                        elseif i1 == 2
                            xx = 2*xu-1;    yy = 2*yu;
                        elseif i1 == 3
                            xx = 2*xu;      yy = 2*yu-1;
                        else
                            xx = 2*xu;      yy = 2*yu;
                        end
                        if BL(xx,yy) ~= 7
                            des = 0;
                        end
                    end
                    if des == 1
                        BL(xu,yu) = 7;
                    end
                    if n > 1
                        UB = UB(1:end-1,1:end);
                    else
                        UB = [];
                        if sub < 3
                            sub = sub + 1;
                        else
                            sub = 1;
                            if j < C/(2^(L+1))
                                j = j + 1;
                            else
                                j = 1;
                                i = i + 1;
                            end
                        end
                    end
                end
            end
            %%%%%%%% Bitstream Stream Truncation %%%%%%%%%%%%
            if bitcount < rate*size(S,1)*size(S,2) || T > 0.5
                T = T/2;
                I = 1;
            else
                I = 50;
            end
            last = find(BS,1,'last');
            if numel(BS) == last
                nz = numel(BS) - last;
                BS = BS(1:last);
            elseif numel(BS) == last + 1
                nz = numel(BS) - last-1;
                BS = BS(1:last+1);
            else
                nz = numel(BS) - last-2;
                BS = BS(1:last+2);
            end
            Stream{iter} = BS;
            iter = iter + 1;
            bitcount = bitcount - nz;
            BS = [];
        else
            break
        end
    end
    if numel(BS) ~= 0
        last = find(BS,1,'last');
        if numel(BS) == last
            nz = numel(BS) - last;
            BS = BS(1:last);
        elseif numel(BS) == last + 1
            nz = numel(BS) - last-1;
            BS = BS(1:last+1);
        else
            nz = numel(BS) - last-2;
            BS = BS(1:last+2);
        end
        Stream{iter} = BS;
        iter = iter + 1;
        bitcount = bitcount - nz;
        BS = [];
    end