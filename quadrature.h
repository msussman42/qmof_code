void init_quadrature( int q )
{
    switch( q )
    {
    case(3):
        qs[0] = 0.1127016653792583;
        qs[1] = 0.5;
        qs[2] = 0.8872983346207417;

        ws[0] = 0.27777777777777785;
        ws[1] = 0.4444444444444444;
        ws[2] = 0.27777777777777785;
        break;
    case(4):
        qs[0] = 0.06943184420297371;
        qs[1] = 0.33000947820757187;
        qs[2] = 0.6699905217924281;
        qs[3] = 0.9305681557970262;

        ws[0] = 0.17392742256872684;
        ws[1] = 0.3260725774312731;
        ws[2] = 0.3260725774312731;
        ws[3] = 0.17392742256872684;
        break;
    case(5):
        qs[0] = 0.04691007703066802;
        qs[1] = 0.23076534494715845;
        qs[2] = 0.5;
        qs[3] = 0.7692346550528415;
        qs[4] = 0.9530899229693319;

        ws[0] = 0.11846344252809471;
        ws[1] = 0.2393143352496831;
        ws[2] = 0.2844444444444445;
        ws[3] = 0.2393143352496831;
        ws[4] = 0.11846344252809471;
        break;
    case(6):
        qs[0] = 0.033765242898423975;
        qs[1] = 0.16939530676686776;
        qs[2] = 0.3806904069584015;
        qs[3] = 0.6193095930415985;
        qs[4] = 0.8306046932331322;
        qs[5] = 0.966234757101576;

        ws[0] = 0.08566224618958487;
        ws[1] = 0.18038078652406947;
        ws[2] = 0.23395696728634569;
        ws[3] = 0.23395696728634569;
        ws[4] = 0.18038078652406947;
        ws[5] = 0.08566224618958487;
        break;
    case(7):
        qs[0] = 0.025446043828620757;
        qs[1] = 0.12923440720030277;
        qs[2] = 0.2970774243113014;
        qs[3] = 0.5;
        qs[4] = 0.7029225756886985;
        qs[5] = 0.8707655927996972;
        qs[6] = 0.9745539561713792;

        ws[0] = 0.06474248308443532;
        ws[1] = 0.1398526957446383;
        ws[2] = 0.19091502525255916;
        ws[3] = 0.20897959183673448;
        ws[4] = 0.19091502525255916;
        ws[5] = 0.1398526957446383;
        ws[6] = 0.06474248308443532;
        break;
    case(8):
        qs[0] = 0.019855071751231912;
        qs[1] = 0.10166676129318664;
        qs[2] = 0.2372337950418355;
        qs[3] = 0.4082826787521751;
        qs[4] = 0.5917173212478248;
        qs[5] = 0.7627662049581645;
        qs[6] = 0.8983332387068134;
        qs[7] = 0.9801449282487681;

        ws[0] = 0.050614268145188344;
        ws[1] = 0.11119051722668717;
        ws[2] = 0.15685332293894352;
        ws[3] = 0.18134189168918088;
        ws[4] = 0.18134189168918088;
        ws[5] = 0.15685332293894352;
        ws[6] = 0.11119051722668717;
        ws[7] = 0.050614268145188344;
        break;
    case(9):
        qs[0] = 0.015919880246186957;
        qs[1] = 0.08198444633668212;
        qs[2] = 0.19331428364970482;
        qs[3] = 0.33787328829809554;
        qs[4] = 0.5;
        qs[5] = 0.6621267117019045;
        qs[6] = 0.8066857163502952;
        qs[7] = 0.9180155536633179;
        qs[8] = 0.984080119753813;

        ws[0] = 0.04063719418078736;
        ws[1] = 0.09032408034742856;
        ws[2] = 0.13030534820146783;
        ws[3] = 0.1561735385200014;
        ws[4] = 0.16511967750062984;
        ws[5] = 0.1561735385200014;
        ws[6] = 0.13030534820146783;
        ws[7] = 0.09032408034742856;
        ws[8] = 0.04063719418078736;
        break;
    case(10):
        qs[0] = 0.013046735741414128;
        qs[1] = 0.06746831665550773;
        qs[2] = 0.16029521585048778;
        qs[3] = 0.2833023029353764;
        qs[4] = 0.4255628305091844;
        qs[5] = 0.5744371694908156;
        qs[6] = 0.7166976970646236;
        qs[7] = 0.8397047841495122;
        qs[8] = 0.9325316833444923;
        qs[9] = 0.9869532642585859;

        ws[0] = 0.033335672154344034;
        ws[1] = 0.07472567457529018;
        ws[2] = 0.109543181257991;
        ws[3] = 0.13463335965499826;
        ws[4] = 0.1477621123573765;
        ws[5] = 0.1477621123573765;
        ws[6] = 0.13463335965499826;
        ws[7] = 0.109543181257991;
        ws[8] = 0.07472567457529018;
        ws[9] = 0.033335672154344034;
        break;
    case(11):
        qs[0] = 0.010885670926971514;
        qs[1] = 0.05646870011595234;
        qs[2] = 0.13492399721297532;
        qs[3] = 0.2404519353965941;
        qs[4] = 0.36522842202382755;
        qs[5] = 0.5;
        qs[6] = 0.6347715779761725;
        qs[7] = 0.7595480646034058;
        qs[8] = 0.8650760027870247;
        qs[9] = 0.9435312998840477;
        qs[10] = 0.9891143290730284;

        ws[0] = 0.027834283558086582;
        ws[1] = 0.06279018473245235;
        ws[2] = 0.09314510546386721;
        ws[3] = 0.11659688229599534;
        ws[4] = 0.13140227225512338;
        ws[5] = 0.13646254338895045;
        ws[6] = 0.13140227225512338;
        ws[7] = 0.11659688229599534;
        ws[8] = 0.09314510546386721;
        ws[9] = 0.06279018473245235;
        ws[10] = 0.027834283558086582;
        break;
    case(12):
        qs[0] = 0.009219682876640378;
        qs[1] = 0.0479413718147626;
        qs[2] = 0.11504866290284765;
        qs[3] = 0.20634102285669126;
        qs[4] = 0.31608425050090994;
        qs[5] = 0.43738329574426554;
        qs[6] = 0.5626167042557344;
        qs[7] = 0.6839157494990901;
        qs[8] = 0.7936589771433087;
        qs[9] = 0.8849513370971523;
        qs[10] = 0.9520586281852375;
        qs[11] = 0.9907803171233596;

        ws[0] = 0.02358766819325601;
        ws[1] = 0.05346966299765944;
        ws[2] = 0.08003916427167306;
        ws[3] = 0.10158371336153282;
        ws[4] = 0.11674626826917732;
        ws[5] = 0.12457352290670134;
        ws[6] = 0.12457352290670134;
        ws[7] = 0.11674626826917732;
        ws[8] = 0.10158371336153282;
        ws[9] = 0.08003916427167306;
        ws[10] = 0.05346966299765944;
        ws[11] = 0.02358766819325601;
        break;
    case(13):
        qs[0] = 0.007908472640705932;
        qs[1] = 0.04120080038851104;
        qs[2] = 0.09921095463334506;
        qs[3] = 0.17882533027982989;
        qs[4] = 0.2757536244817766;
        qs[5] = 0.3847708420224326;
        qs[6] = 0.5;
        qs[7] = 0.6152291579775674;
        qs[8] = 0.7242463755182234;
        qs[9] = 0.8211746697201701;
        qs[10] = 0.9007890453666549;
        qs[11] = 0.958799199611489;
        qs[12] = 0.9920915273592941;

        ws[0] = 0.02024200238265794;
        ws[1] = 0.0460607499188643;
        ws[2] = 0.06943675510989368;
        ws[3] = 0.08907299038097276;
        ws[4] = 0.10390802376844428;
        ws[5] = 0.11314159013144857;
        ws[6] = 0.11627577661543695;
        ws[7] = 0.11314159013144857;
        ws[8] = 0.10390802376844428;
        ws[9] = 0.08907299038097276;
        ws[10] = 0.06943675510989368;
        ws[11] = 0.0460607499188643;
        ws[12] = 0.02024200238265794;
        break;
    case(14):
        qs[0] = 0.006858095651593843;
        qs[1] = 0.03578255816821324;
        qs[2] = 0.08639934246511749;
        qs[3] = 0.15635354759415726;
        qs[4] = 0.24237568182092295;
        qs[5] = 0.3404438155360551;
        qs[6] = 0.44597252564632817;
        qs[7] = 0.5540274743536718;
        qs[8] = 0.6595561844639448;
        qs[9] = 0.757624318179077;
        qs[10] = 0.8436464524058427;
        qs[11] = 0.9136006575348825;
        qs[12] = 0.9642174418317868;
        qs[13] = 0.9931419043484062;

        ws[0] = 0.017559730165876187;
        ws[1] = 0.04007904357988015;
        ws[2] = 0.06075928534395148;
        ws[3] = 0.0786015835790967;
        ws[4] = 0.09276919873896881;
        ws[5] = 0.10259923186064777;
        ws[6] = 0.10763192673157883;
        ws[7] = 0.10763192673157883;
        ws[8] = 0.10259923186064777;
        ws[9] = 0.09276919873896881;
        ws[10] = 0.0786015835790967;
        ws[11] = 0.06075928534395148;
        ws[12] = 0.04007904357988015;
        ws[13] = 0.017559730165876187;
        break;
    case(15):
        qs[0] = 0.006003740989757311;
        qs[1] = 0.031363303799647024;
        qs[2] = 0.0758967082947864;
        qs[3] = 0.13779113431991497;
        qs[4] = 0.21451391369573058;
        qs[5] = 0.3029243264612183;
        qs[6] = 0.39940295300128276;
        qs[7] = 0.5;
        qs[8] = 0.6005970469987173;
        qs[9] = 0.6970756735387817;
        qs[10] = 0.7854860863042694;
        qs[11] = 0.862208865680085;
        qs[12] = 0.9241032917052137;
        qs[13] = 0.968636696200353;
        qs[14] = 0.9939962590102427;

        ws[0] = 0.015376620998059323;
        ws[1] = 0.035183023744054034;
        ws[2] = 0.053579610233585886;
        ws[3] = 0.06978533896307695;
        ws[4] = 0.08313460290849689;
        ws[5] = 0.09308050000778094;
        ws[6] = 0.09921574266355562;
        ws[7] = 0.10128912096278045;
        ws[8] = 0.09921574266355562;
        ws[9] = 0.09308050000778094;
        ws[10] = 0.08313460290849689;
        ws[11] = 0.06978533896307695;
        ws[12] = 0.053579610233585886;
        ws[13] = 0.035183023744054034;
        ws[14] = 0.015376620998059323;
        break;
    case(16):
        qs[0] = 0.005299532504175031;
        qs[1] = 0.0277124884633837;
        qs[2] = 0.06718439880608412;
        qs[3] = 0.1222977958224985;
        qs[4] = 0.19106187779867811;
        qs[5] = 0.2709916111713863;
        qs[6] = 0.35919822461037054;
        qs[7] = 0.4524937450811813;
        qs[8] = 0.5475062549188188;
        qs[9] = 0.6408017753896295;
        qs[10] = 0.7290083888286136;
        qs[11] = 0.8089381222013219;
        qs[12] = 0.8777022041775016;
        qs[13] = 0.9328156011939159;
        qs[14] = 0.9722875115366163;
        qs[15] = 0.994700467495825;

        ws[0] = 0.013576229705877019;
        ws[1] = 0.031126761969323853;
        ws[2] = 0.047579255841246296;
        ws[3] = 0.062314485627767015;
        ws[4] = 0.07479799440828838;
        ws[5] = 0.08457825969750131;
        ws[6] = 0.0913017075224618;
        ws[7] = 0.09472530522753429;
        ws[8] = 0.09472530522753429;
        ws[9] = 0.0913017075224618;
        ws[10] = 0.08457825969750131;
        ws[11] = 0.07479799440828838;
        ws[12] = 0.062314485627767015;
        ws[13] = 0.047579255841246296;
        ws[14] = 0.031126761969323853;
        ws[15] = 0.013576229705877019;
        break;
    default:
        printf("DONT HAVE CORRECT QUADRATURE\n");
    }
}