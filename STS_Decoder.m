%% STS Decoder %% 
function [Rec_S,n_sig_coeff] = STS_Decoder(Stream,row,col,L,Tint)
Rec_S = zeros(row,col);tic
[R,C] = size(Rec_S);
BL = zeros(R/2,C/2);
[br,bc] = size(BL);
BL(br/2+1:br,1:bc) = 3;
BL(1:br/2,bc/2+1:bc) = 3;
BL(1:R/(2^(L)),1:C/(2^(L))) = 4;
T = Tint;mx = 1;a1 = (R/(2^(L+1))); l1 = 0;iter = 1;
% Entropy Decoder
while iter <= size(Stream,2)
Str = [];
Str = Stream{mx};
Nstr = numel(Str);n1 = find(Str == 1, 1);
sg = 1;m1 = 0;a = 0;o1 = 0;i = 1; % from here, i tracks bits in str
[x,y] = find(BL > 3);
    while sg <= 4*numel(x) % sg gives the pixel number in refinement pass
        if i > Nstr
            break
        end
        n = ceil(sg/4); % x(n),y(n) gives BL location
        r = rem(sg,4);
        if r == 0
            r = 4;
        end; 
        xx = ceil(r/2); % gives x coordinate location within BL location
        yy = rem(r,2); % gives y coordinate location within BL location
        if yy == 0
            yy = 2;
        end
         p = 2*x(n)-2+xx; q = 2*y(n)-2+yy;
         if Rec_S(p,q) == 0
             if sg < n1
                 sg = sg + 1;
                 i = i + 1;
             else
                if Str(i) == 0 && a == 0
                    sg = sg + 1;
                    i = i + 1;
                elseif Str(i) == 0 && a == 1
                    sg = sg + 1;
                    i = i + 1;
                    Rec_S(p,q) = 1.5*T;
                    a = 0;
                elseif Str(i) == 1 && a == 0
                    a = 1;
                    i = i + 1;
                elseif Str(i) == 1 && a == 1
                    sg = sg+1;
                    i = i + 1;
                    Rec_S(p,q) = - 1.5*T;
                    a = 0;
                end;
             end
         elseif Rec_S(p,q) < 0
             if sg < n1
                 sg = sg + 1;
                 i = i + 1;
                 Rec_S(p,q) = Rec_S(p,q) + 0.5*T;
             else
                if Str(i) == 1
                    Rec_S(p,q) = Rec_S(p,q) - 0.5*T;
                    sg = sg + 1;
                    i = i + 1;
                else
                    Rec_S(p,q) = Rec_S(p,q) + 0.5*T;
                    sg = sg + 1;
                    i = i + 1;
                end
             end
         elseif Rec_S(p,q) > 0
             if sg < n1
                 sg = sg + 1;
                 i = i + 1;
                 Rec_S(p,q) = Rec_S(p,q) - 0.5*T;
             else
                if Str(i) == 1
                    Rec_S(p,q) = Rec_S(p,q) + 0.5*T;
                    sg = sg + 1;
                    i = i + 1;
                else
                    Rec_S(p,q) = Rec_S(p,q) - 0.5*T;
                    sg = sg + 1;
                    i = i + 1;
                end
             end
         end
    end
    p = 1;
 while p <= R/(2^(L+1)) % SIGNIFICANCE PASS
     if i > Nstr
         break
     end
     q = 1;sub = 1;
     while q <= C/(2^(L+1))
         if i > Nstr
             break
         end
         if sub == 1 %LH
             BL_addr_x = p; BL_addr_y = q + C/(2^(L+1)); % LH
         elseif sub == 2 %HL
             BL_addr_x = p + R/(2^(L+1)); BL_addr_y = q; % HL
         else %HH
             BL_addr_x = p + R/(2^(L+1)); BL_addr_y = q + C/(2^(L+1)); % HH
         end
         if BL(BL_addr_x, BL_addr_y) == 7
             if sub < 3
                 sub = sub + 1;
             else
                 sub = 1;
                 q = q + 1;
             end
         elseif BL(BL_addr_x,BL_addr_y) == 4
             UB = [BL_addr_x,BL_addr_y];
             LB = [2*BL_addr_x-1, 2*BL_addr_y-1;2*BL_addr_x-1, 2*BL_addr_y;...
                 2*BL_addr_x, 2*BL_addr_y-1;2*BL_addr_x, 2*BL_addr_y];
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
         sz = size(LB,1);l1 = 0;
         while sz ~= 0
             if i > Nstr
                 break
             end
             if BL(LB(1,1),LB(1,2)) <= 3
                 % if block insignificant
                 if l1 == 0
                     if Str(i) == 0
                        if sz > 1
                            LB = LB(2:end,:);
                        else
                            LB = [];
                        end
                        i = i + 1;
                        sz = sz - 1;
                     else
                        i = i + 1;
                        l1 = 1;
                        m1 = 0; o1 = 0;
                     end
                 elseif l1 == 1;
                     if m1 == 0 % m1 tracks the actual coefficient number in BL
                        if Str(i) == 0 && o1 == 0
                            m1 = 1;
                            if i == Nstr
                                BL(LB(1,1),LB(1,2)) = BL(LB(1,1),LB(1,2)) + 4;
                            end
                            i = i + 1;
                        elseif Str(i) == 0 && o1 == 1
                            Rec_S(2*LB(1,1)-1,2*LB(1,2)-1) = 1.5*T;
                            m1 = 1;
                            if i == Nstr
                                BL(LB(1,1),LB(1,2)) = BL(LB(1,1),LB(1,2)) + 4;
                            end
                            i = i + 1;
                            o1 = 0;
                        elseif Str(i) == 1 && o1 == 0
                            o1 = 1;
                            i = i + 1;
                        elseif Str(i) == 1 && o1 == 1
                            Rec_S(2*LB(1,1)-1,2*LB(1,2)-1) = -1.5*T;
                            m1 = 1;
                            if i == Nstr
                                BL(LB(1,1),LB(1,2)) = BL(LB(1,1),LB(1,2)) + 4;
                            end
                            i = i + 1;
                            o1 = 0;
                        end
                    elseif m1 == 1 % m tracks the actual coefficient number in BL
                        if Str(i) == 0 && o1 == 0
                            m1 = 2;
                            if i == Nstr
                                BL(LB(1,1),LB(1,2)) = BL(LB(1,1),LB(1,2)) + 4;
                            end
                            i = i + 1;
                        elseif Str(i) == 0 && o1 == 1
                            Rec_S(2*LB(1,1)-1,2*LB(1,2)) = 1.5*T;
                            m1 = 2;
                            if i == Nstr
                                BL(LB(1,1),LB(1,2)) = BL(LB(1,1),LB(1,2)) + 4;
                            end
                            i = i + 1;
                            o1 = 0;
                        elseif Str(i) == 1 && o1 == 0
                            o1 = 1;
                            i = i + 1;
                        elseif Str(i) == 1 && o1 == 1
                            Rec_S(2*LB(1,1)-1,2*LB(1,2)) = -1.5*T;
                            m1 = 2;
                            if i == Nstr
                                BL(LB(1,1),LB(1,2)) = BL(LB(1,1),LB(1,2)) + 4;
                            end
                            i = i + 1;
                            o1 = 0;
                        end
                    elseif m1 == 2 % m tracks the actual coefficient number in BL
                        if Str(i) == 0 && o1 == 0
                            m1 = 3;
                            if i == Nstr
                                BL(LB(1,1),LB(1,2)) = BL(LB(1,1),LB(1,2)) + 4;
                            end
                            i = i + 1;
                        elseif Str(i) == 0 && o1 == 1
                            Rec_S(2*LB(1,1),2*LB(1,2)-1) = 1.5*T;
                            m1 = 3;
                            if i == Nstr
                                BL(LB(1,1),LB(1,2)) = BL(LB(1,1),LB(1,2)) + 4;
                            end
                            i = i + 1;
                            o1 = 0;
                        elseif Str(i) == 1 && o1 == 0
                            o1 = 1;
                            i = i + 1;
                        elseif Str(i) == 1 && o1 == 1
                            Rec_S(2*LB(1,1),2*LB(1,2)-1) = -1.5*T;
                            m1 = 3;
                            if i == Nstr
                                BL(LB(1,1),LB(1,2)) = BL(LB(1,1),LB(1,2)) + 4;
                            end
                            i = i + 1;
                            o1 = 0;
                        end
                    elseif m1 == 3 % m tracks the actual coefficient number in BL
                        if Str(i) == 0 && o1 == 0
                            m1 = 0;l1 = 0;
                            i = i + 1;
                            BL(LB(1,1),LB(1,2)) = BL(LB(1,1),LB(1,2)) + 4;
                        elseif Str(i) == 0 && o1 == 1
                            Rec_S(2*LB(1,1),2*LB(1,2)) = 1.5*T;
                            m1 = 0;
                            i = i + 1;
                            o1 = 0;l1 = 0;
                            BL(LB(1,1),LB(1,2)) = BL(LB(1,1),LB(1,2)) + 4;
                        elseif Str(i) == 1 && o1 == 0
                            o1 = 1;
                            i = i + 1;
                        elseif Str(i) == 1 && o1 == 1
                            Rec_S(2*LB(1,1),2*LB(1,2)) = -1.5*T;
                            m1 = 0;
                            i = i + 1;
                            o1 = 0;l1 = 0;
                            BL(LB(1,1),LB(1,2)) = BL(LB(1,1),LB(1,2)) + 4;
                        end
                     end
                 end
                 
             else
                 % current block significant, now check for descendants
                 if BL(LB(1,1),LB(1,2)) == 4
                     if Str(i) == 1
                         i = i + 1;
                         x1 = 2*LB(1,1)-1;y1 = 2*LB(1,2)-1;
                         x2 = 2*LB(1,1)-1;y2 = 2*LB(1,2);
                         x3 = 2*LB(1,1);y3 = 2*LB(1,2)-1;
                         x4 = 2*LB(1,1);y4 = 2*LB(1,2);
                         BL(LB(1,1),LB(1,2)) = BL(LB(1,1),LB(1,2)) + 1;
                         UB = [UB; LB(1,1),LB(1,2)];
                         if sz > 1
                            LB = LB(2:end,:);
                         else
                            LB = [];
                         end
                         LB = [x1, y1; x2, y2; x3, y3; x4, y4; LB];
                         sz = sz + 3;
                     else
                         i = i + 1;
                         if sz > 1
                            LB = LB(2:end,1:end);
                         else
                            LB = [];
                         end
                         sz = sz - 1;
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
                             LB1 = [LB1; x1, y1];
                             sz = sz + 1;
                         end
                     end
                     UB = [UB; LB(1,1),LB(1,2)];
                     if sz > 1
                         LB = LB(2:end,1:end);
                     else
                         LB = [];
                     end
                     sz = sz - 1;
                     LB = [LB1; LB];
                 else
                     if sz > 1
                         LB = LB(2:end,1:end);
                     else
                         LB = [];
                     end
                     sz = sz - 1;
                 end
             end
         end
         
         % UB Part LIFO
         szu = size(UB,1);
         while szu ~= 0
             xu = UB(szu,1);   yu = UB(szu,2);
             des = 1; tpb = 0;
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
                 else
                     tpb = 1;
                 end
             end
             if BL(xu,yu) == 4
                 if des == 1
                    BL(xu,yu) = 7;
                 elseif tpb == 0
                     BL(xu,yu) = 4;
                 else
                     BL(xu,yu) = 5;
                 end
             elseif BL(xu,yu) == 5
                 if des == 1
                     BL(xu,yu) = 7;
                 end
             end
             if size(UB,1) > 1
                 UB = UB(1:end-1,1:end);
             else
                 UB = [];
             end
             szu = szu -1;
         end
         if sub < 3
             sub = sub + 1;
         else
             sub = 1;
             q = q + 1;
         end
     end
     p = p + 1;
 end
% Parameter Evaluation Part
iter = iter + 1;
T = T/2;
mx = mx + 1;
n_sig_coeff(iter-1) = numel(find(Rec_S ~= 0));
end
