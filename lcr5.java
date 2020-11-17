/**
 *
 * @author rowew
 */
import java.math.*;
import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.Date;

public class lcr5 {
    BufferedReader inputReader =null;
    BufferedWriter outputWriter =null;
    public lcr5()
           
    {
        SimpleDateFormat formatter= new SimpleDateFormat("yy-MM-dd");
Date date = new Date(System.currentTimeMillis());
String holder = String.valueOf(formatter.format(date));
//System.out.println(holder);
String splitter [] = holder.split("-");
String fileName =splitter[0]+splitter[1]+splitter[2]+"LCR";
//System.out.println(fileName);
     
        Vector lines = new Vector();
        int total1=0;
        int total2 =0;
        String location ="out//"+fileName+"//6//plates//";
        String location2 ="out//"+fileName+"//6//";
       String location3 ="out//"+fileName+"//4//";
        String lett ="ABCDEFGHIJKLMNOPQRSTUVW";
        String plate[][]= new String[19][19];
        Vector temp = new Vector();
        for(int i=0; i< plate.length; i++)
        {
            for(int j=0; j< plate.length; j ++)
            {
                plate[i][j] ="";
            }
        }
     for(int i=1; i< plate.length; i++)
     {
       plate[0][i] ="id";
       plate[1][i] =String.valueOf(i);
       if(i >1)
       {
           plate[i][0]= String.valueOf(lett.charAt(i-2));
       }
         
     }
//System.out.println(location);
 File folder = new File(location);
File[] listOfFiles = folder.listFiles();
String name="";
//System.out.println(location+name);


for (File file : listOfFiles) {
    if (file.isFile() && String.valueOf(file).indexOf("POOLCR") >-1) {
     name= file.getName();
       
    }
}
        int count=0;
       fileReader(location+name);
        String currentLine =null;
       
           while(true)
        {
          try{

                 currentLine = inputReader.readLine();
          }
          catch(Exception e)
          {
          }
          if(currentLine==null)
          {
              break;
          }
          count++;
         if(count <3)
         {
           
             continue;
         }
       
          String splits [] =currentLine.split(",");
          for(int i=1; i< splits.length; i++)
          {
         plate[count-1][i] =splits[i];
         if(splits[i].indexOf("_dig") >0)
         {
               if(i >total1)
          {
              total1 =i;
          }
            if(count-1 >total2)
            {
                total2 =count-1;
            }
         }
                 }
       
         
        }
           
           
for (File file : listOfFiles) {
    if (file.isFile() && String.valueOf(file).indexOf("LCRLCR") >-1) {
     name= file.getName();
       
    }
}
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
          count=0;
       fileReader(location+name);
         currentLine =null;
    //    System.out.println(total1+" "+total2);
           while(true)
        {
          try{

                 currentLine = inputReader.readLine();
          }
          catch(Exception e)
          {
          }
          if(currentLine==null)
          {
              break;
          }
          count++;
         if(count <3)
         {
           
             continue;
         }
       
          String splits [] =currentLine.split(",");
          for(int i=1; i< splits.length; i++)
          {
             
              if(splits[i].startsWith("SBC"))
              {
                //  System.out.println(count);
         plate[count-1][i+total1+2] =splits[i];
              }
                 }
       
         
        }


           plate[total2+2][1]="water";
           plate[total2+3][1]="ampligase";
           plate[total2+4][1]="mm_lcr";          
           
           fileWriter("out//"+fileName+"//"+String.valueOf(name.substring(0,6))+"_Echo_source_plate.csv");
         
           for(int i=0; i< plate.length; i++)
           {
               for(int j=0; j< plate.length; j++)
               {
                //  System.out.println(plate[i][j]);
                   sequenceWriter(plate[i][j]);
                   
                      if(j <plate.length-1)
                      {
                          sequenceWriter(",");
                      }
               }
               sequenceWriter("\n");
               
           }
           try{
               outputWriter.close();
           }
           catch(Exception e)
           {
               
           }
           
           
   
   
   
   
   
         count=0;
       fileReader(location+name);
         currentLine =null;
    //    System.out.println(total1+" "+total2);
           while(true)
        {
          try{

                 currentLine = inputReader.readLine();
          }
          catch(Exception e)
          {
          }
          if(currentLine==null)
          {
              break;
          }
          count++;
         if(count <3)
         {
           
             continue;
         }
       
          String splits [] =currentLine.split(",");
          for(int i=1; i< splits.length; i++)
          {
             
              if(splits[i].startsWith("SBC"))
              {
                //  System.out.println(count);
         plate[count-1][i+total1+2] =splits[i];
              }
                 }
       
         
        }
           
           
           double totWat =0.0;
double totMM =0.0;
double totAmp =0.0;



int y =total2+5;
int x =1;

 count =0;
         
           fileReader("out//"+fileName+"//6//worklist.csv");
              while(true)
        {
          try{

                 currentLine = inputReader.readLine();
          }
          catch(Exception e)
          {
          }
          if(currentLine==null)
          {
              break;
          }
         count++;
         if(count ==1)
         {
         //    lines.add(new String(currentLine));
             continue;
         }
               
               String splits [] =currentLine.split(",");
               
               if(splits[5].equals("water"))
               {
                 totWat +=Double.parseDouble(splits[0]);
                 continue;
               }
                  if(splits[5].equals("mm_lcr"))
               {
                 totMM +=Double.parseDouble(splits[0]);
               continue;
               }
               
                     if(splits[5].equals("ampligase"))
               {
                 totAmp +=Double.parseDouble(splits[0]);
                 continue;
               }
               
               if(splits[5].indexOf("_dominoes") > -1)
               {
                   plate[y][x] =splits[5];
                   x++;
                   if(x > 18)
                   {
                       x=1;
                       y++;
                   }
                   
               }
                     
                     
             
           String out2="";
           
       
        }
           


//System.out.println(totMM);
//System.out.println(totWat);
//System.out.println(totAmp);


int mmTot =(int)(totMM/50.0)+2;
int watTot =(int)(totWat/50.0)+2;
int ampTot =(int)(totAmp/50.0)+2;

for(int i=0; i< watTot; i++)
{
           plate[total2+2][i+1]="water";
}
for(int i=0; i < ampTot; i++)
{
           plate[total2+3][1+i]="ampligase";
}
for(int i=0; i< mmTot; i++)
{
           plate[total2+4][1+i]="mm_lcr";    
}
           
           fileWriter("out//"+fileName+"//"+String.valueOf(name.substring(0,6))+"_Echo_source_plate.csv");
         
           for(int i=0; i< plate.length; i++)
           {
               for(int j=0; j< plate.length; j++)
               {
                //  System.out.println(plate[i][j]);
                   sequenceWriter(plate[i][j]);
                   
                      if(j <plate.length-1)
                      {
                          sequenceWriter(",");
                      }
               }
               sequenceWriter("\n");
               
           }
           try{
               outputWriter.close();
           }
           catch(Exception e)
           {
               
           }



















           
           int cWat=1;
           int cAmp=1;
           int cMM =1;
           
           double curWat=0.0;
           double curAmp =0.0;
           double curMM =0.0;






            count =0;
           fileWriter("out//"+fileName+"//"+String.valueOf(name.substring(0,6))+"_Echo_worklist.csv");
           sequenceWriter("Source Plate Name,Source Plate Type,Component,Source Well,Destination Plate Name,Destination Well,Transfer Volume,Destination Well X Offset,Destination Well Y Offset,Delay,Destination name\n");
           
           fileReader("out//"+fileName+"//6//worklist.csv");
              while(true)
        {
          try{

                 currentLine = inputReader.readLine();
          }
          catch(Exception e)
          {
          }
          if(currentLine==null)
          {
              break;
          }
         count++;
         if(count ==1)
         {
             lines.add(new String(currentLine));
             continue;
         }
               
               String splits [] =currentLine.split(",");
               
               String out=String.valueOf(name.substring(0,6))+" Echo source,384PP_AQ_GP3,"+splits[5]+",";
               String source =getSource(splits[5], plate);
             
               
            //   System.out.println(splits[5]+" "+source);
               out+=source+",";
               out+=String.valueOf(name.substring(0,6))+"ECHLCR,";
               out+=splits[4]+",";
               double hold = Double.parseDouble(splits[0])*100;
               out+=hold+",,,,";
               out+=splits[6];
        //System.out.println(out);
        splits[4]= source;
        if(currentLine.indexOf("water") >-1)
        {
            double cur = Double.parseDouble(splits[0]);
            if(curWat+cur >50)
            {
                cWat++;
                curWat=cur;
            }
            else{
                curWat+=cur;
            }
            splits[4] =String.valueOf(source.charAt(0))+cWat;
        }
        if(currentLine.indexOf("mm_lcr") >-1)
        {
             double cur = Double.parseDouble(splits[0]);
            if(curMM+cur >50)
            {
                cMM++;
                curMM=cur;
            }
            else{
                curMM+=cur;
            }
            splits[4] =String.valueOf(source.charAt(0))+cMM;
        }
        if(currentLine.indexOf("ampligase") >-1)
        {
                double cur = Double.parseDouble(splits[0]);
            if(curAmp+cur >50)
            {
                cAmp++;
                curAmp=cur;
            }
            else{
                curAmp+=cur;
            }
            splits[4] =String.valueOf(source.charAt(0))+cAmp;
        }
        splits[6] =splits[5];
           sequenceWriter(out+"\n");
           String out2="";
           for(int i=0 ;i < splits.length; i++)
           {
              out2+=splits[i]+",";
           }
           out2 =out2.substring(0, out2.length()-1);
               lines.add(new String(out2));
        }
             
             
             
             
             
             
              try{
                 
                  outputWriter.close();
              }
              catch(Exception e)
              {
                 
              }
             
               fileWriter("out//"+fileName+"//"+String.valueOf(name.substring(0,6))+"_Echo_source_plate_worklist.csv");
     
            int countAmp=0;
            int countWat=0;
            int countMM=0;
            sequenceWriter(String.valueOf(lines.get(0))+"\n");
              for(int i=1; i< lines.size(); i++)
              {
                  String line =String.valueOf(lines.get(i));
           
         
               
                  if(line.indexOf("water") >-1)
                          {
                         
                        if(countWat ==0)
                        {
                          //  System.out.println("yesss");
                            countWat++;
                           
                            for(int k=0; k < watTot; k++)
                            {
                                      String splits [] =line.split(",");
                                 splits[0] ="65";
                                 String out2 ="";
                                 splits[4] =String.valueOf(splits[4].charAt(0))+(k+1);
                   for(int j=0; j< splits.length; j++)
                   {
                      out2+=splits[j]+",";
                   }
                   out2=out2.substring(0, out2.length()-1);
                   sequenceWriter(out2+"\n");
                            }
                           
                        }
                   
                           }
                         
               else     if(line.indexOf("mm_lcr") >-1)
                          {
                             
                        if(countMM ==0)
                        {
                            countMM++;
                           
                            for(int k=0; k < mmTot; k++)
                            {
                                      String splits [] =line.split(",");
                                 splits[0] ="65";
                                 String out2 ="";
                                // System.out.println(line);
                              //  System.out.println(splits[4]);
                                 splits[4] =String.valueOf(splits[4].charAt(0))+(k+1);
                   for(int j=0; j< splits.length; j++)
                   {
                      out2+=splits[j]+",";
                   }
                   out2=out2.substring(0, out2.length()-1);
                   sequenceWriter(out2+"\n");
                            }
                           
                        }
                          }
                     else if(line.indexOf("ampligase") >-1)
                          {
                                         if(countAmp ==0)
                        {
                            countAmp++;
                           
                            for(int k=0; k < ampTot; k++)
                            {
                                      String splits [] =line.split(",");
                                 splits[0] ="65";
                                 String out2 ="";
                                 splits[4] =String.valueOf(splits[4].charAt(0))+(k+1);
                   for(int j=0; j< splits.length; j++)
                   {
                      out2+=splits[j]+",";
                   }
                   out2=out2.substring(0, out2.length()-1);
                   sequenceWriter(out2+"\n");
                            }
                           
                        }
                          }
                     
                     else{
                         
                    String splits[] = line.split(",");
                    String out2="";
                    splits[0] ="65";

                   for(int j=0; j< splits.length; j++)
                   {
                      out2+=splits[j]+",";
                   }
                   out2=out2.substring(0, out2.length()-1);
                   sequenceWriter(out2+"\n");

                //  sequenceWriter(String.valueOf(lines.get(i))+"\n");
                     }
              }
              try{
                  outputWriter.close();
              }
              catch(Exception e)
              {
                 
              }
           
    }
    public static void main (String args[])
    {
        lcr5 lcr5 = new lcr5 ();
    }
        public void fileReader(String fileName)
    {
         //    System.out.println(fileName);
    try
      { inputReader =
        new BufferedReader
        ( new InputStreamReader
          ( new FileInputStream (fileName)));
      }
      catch (FileNotFoundException f)
        { System.out.println(fileName+" File not found");
        System.exit(1);
        }

    }
        public void sequenceWriter(String lines)
     {
         try{
             outputWriter.write(lines);
        //     outputWriter.newLine();

         }
         catch(Exception e)
         {
             System.out.println(e);
             System.exit(1);
         }
     }

                                             public void fileWriter(String called)
  {


  File file1 = new File(called);
    try
    {
     outputWriter = new BufferedWriter ( new OutputStreamWriter ( new FileOutputStream (file1 )));
    }
    catch (FileNotFoundException f)
    {
      System.out.println("File not created");

      System.exit(2);
    }}
        public String getSource(String source, String[][]plates )
        {
            if(source.indexOf("dominoes") >-1)
            {
             //   source = source.substring(0,source.indexOf("dominoes") -1 );
            }
            String out="";
            for(int i=2; i< plates.length; i++)
            {
                for(int j=1; j < plates.length; j++)
                {
                 if(source.equals(plates[i][j]))
                 {
                  //   System.out.println(plates[j][0]);
                     out=plates[i][0]+plates[1][j];
                 }
                }
            }
         //   System.out.println(out);
            return out;
        }
}

