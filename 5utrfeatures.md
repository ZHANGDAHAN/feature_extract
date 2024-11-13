# feature_extract
ff=fopen('rice_5utr.fasta.out');data=cell(1,1);t=0;
while ~feof(ff)
    f=fgetl(ff);
    t=t+1;data{t,1}=f(2:end);
    f=fgetl(ff);
    data{t,2}=f;
    f=fgetl(ff);
    data{t,3}=f;
    data{t,4}=(f((length(data{t,2})+3):(length(f)-1)));
    data{t,5}=data{t,3}(1:length(data{t,2}));
    f=fgetl(ff);
    f=fgetl(ff);
    data{t,6}=f(1:length(data{t,2})); %% 第六列是high probability structure
    f=fgetl(ff);
end
fclose(ff)
%% 提取特征
TE_mus=zeros(1,1);
tezheng=zeros(1,8);
[uniquedids,a1,b1]=intersect(transcript_mrna(:,1),rawdata(:,1));
transcript_mrnahasutr=transcript_mrna(a1,:);
rawdatahasutr=rawdata(b1,:);
% 数据对齐
[uniquedids,a2,b2]=intersect(rawdatahasutr(:,1),data(:,1));
rawdatahasutr=rawdatahasutr(a2,:);
transcript_mrnahasutr=transcript_mrnahasutr(a2,:);
data=data(b2,:);
for i=1:size(rawdatahasutr,1)
    s=upper(data{i,2});
    % stem区域的序列
    sstem=s([strfind(data{i,5},'('),strfind(data{i,5},')')]);
    % loop区域的序列
    sloop=s(strfind(data{i,5},'.'));
    % stem比例，
    tezheng(i,1)=1-length(strfind(data{i,5},'.'))/length(data{i,5});
    % stem长度，
    tezheng(i,2)=length(data{i,5})-length(strfind(data{i,5},'.'));
    % UTR的长度，
    tezheng(i,3)=length(data{i,5});
    % 自由能，
    tezheng(i,4)=str2double(data{i,4});
    % stem的GC含量，
tezheng(i,5)=(length(strfind(sstem,'G'))+length(strfind(sstem,'C')))/length(sstem);
    % loop的GC含量，
tezheng(i,6)=(length(strfind(sloop,'G'))+length(strfind(sloop,'C')))/length(sloop);
    % 总体的GC含量
    tezheng(i,7)=(length(strfind(s,'G'))+length(strfind(s,'C')))/length(s);
    % UAPE
     if length(data{i,2})>15
        tezheng(i,8)= APE_value(data{i,2}(end-14:end),(rice_UAPE(end-14:end,:)));
    else
        tezheng(i,8)= nan;
     end
    % uORF number download elswhere
end
