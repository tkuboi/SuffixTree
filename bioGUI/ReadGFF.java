package bioGUI;

import java.io.*;
import java.lang.Math;
import java.util.regex.*;
import java.util.Scanner;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Collections;

public class ReadGFF {
    private static ArrayList<mRNA> m_rna;
    private static String m_dna;

    public static ArrayList<mRNA> readfile(String filename) throws IOException {
        m_rna = new ArrayList<mRNA>();
        mRNA rna;
        CDS cds;
        int lines = 0;
        String line = "";
        String geneId = "";
        String transId = "";
        String tokens[];
        boolean direction;
        Scanner s = null;
        System.out.println(filename);
        Pattern p = Pattern.compile("\"(\\w+)\";\\stranscript_id\\s\"(\\S+)\";$");
        Matcher m;

        try {
            s = new Scanner(new BufferedReader(new FileReader(filename)));
            while (s.hasNextLine()) {
                line = s.nextLine();
                if (line.length() <= 0)
                    continue;
                tokens = line.split("\t");
                //for (int i=0;i<tokens.length;i++)
                //    System.out.println(tokens[i]);
                if (tokens[2].equals("mRNA")) {// mRNA
                    rna = new mRNA();
                    rna.setChromosome(tokens[0]);
                    rna.setStart(Integer.parseInt(tokens[3]));
                    rna.setEnd(Integer.parseInt(tokens[4]));
                    rna.setDirection((tokens[6].equals("+")));
                    geneId = "";
                    transId = "";
                    if (tokens.length > 8) {
                        
                        m = p.matcher(tokens[8]);
                        if (m.find()) {
                            geneId = m.group(1);
                            transId = m.group(2);
                        }
                        rna.setGeneID(geneId);
                        rna.setTranscriptID(transId);
                    }
                    //System.out.println(tokens[2]);
                    m_rna.add(rna);
                }
                else {
                    if (tokens[2].equals("CDS")) {//CDS
                        //System.out.println(tokens[2]);
                        m_rna.get(m_rna.size()-1).addCDS(new CDS(Integer.parseInt(tokens[3]),Integer.parseInt(tokens[4]), tokens[6].equals("+")));
                    }
                }
                lines++;
                //System.out.println("---------------------------------");
            }
            if (line.length() == 0)
                lines--;
        } finally {
            if (s != null) {
                s.close();
                System.out.println("lines: " + lines);
            }
        }
        
        Collections.sort(m_rna, new Comparator<mRNA>() {
              public int compare(mRNA mRNAOne, mRNA mRNATwo) {
                 return mRNAOne.getStart() - mRNATwo.getStart();
              }
           }
        );
        
        return m_rna;
    }
    
    public static String readFasta(String filename) throws IOException {
        Scanner s = null;
        int lines = 0;
        String line = "";
        m_dna = "";
        mRNA rna;
        int start, end, idx;
        Pattern p = Pattern.compile("([ATCGN]*)");
        Matcher m;
        ArrayList<CDS> cds;

        try {
            s = new Scanner(new BufferedReader(new FileReader(filename)));
            while (s.hasNextLine()) {
                line = s.nextLine();
                m = p.matcher(line);
                if (m.find())
                    line = m.group(1);
                //System.out.println(line);
                m_dna += line;
            }

        } finally {
            if (s != null) {
                s.close();
            }
        }
        
        /*for (int i = 0; i < m_rna.size(); i++) {
            cds = m_rna.get(i).getCDS();
            end = m_rna.get(i).getEnd();
            start = m_rna.get(i).getStart();
            if (cds.size() > 0) {
                for (int j = 0; j < cds.size(); j++) {
                    end = Math.max(end,cds.get(j).getEnd());
                    start = Math.min(start,cds.get(j).getStart());
                }
            }
            start = start - 1;
            //end = end + 1;
            m_rna.get(i).setSequence(m_dna.substring(start,end));
        }*/
        //System.out.println(m_dna.substring(0,51));
        //return m_rna;
        return m_dna;
    }

    public static ArrayList<mRNA> getmRNA() {
        return m_rna;
    }

    public static void main(String[] args) {
       System.out.println(args[0]);
       System.out.println(args[1]);
       try {
           ReadGFF.readfile(args[0]);
           ReadGFF.readFasta(args[1]);
           ArrayList<mRNA> mRNAs = ReadGFF.getmRNA();
           for (int i = 0; i < mRNAs.size(); i++) {
               System.out.println(mRNAs.get(i).getChromosome() + " " + mRNAs.get(i).getStart() + " " + mRNAs.get(i).getEnd()
                  + " " + mRNAs.get(i).getDirection() + " " + mRNAs.get(i).getGeneID() + " " + mRNAs.get(i).getTranscriptID());
               System.out.println(mRNAs.get(i).getRawSequence());
               System.out.println(mRNAs.get(i).getCorrectedSequence());
               System.out.println(mRNAs.get(i).getCDS().size());
           }
       }
       catch(IOException e) {
           System.out.println("Exception: " + e);
       }
    }
}

