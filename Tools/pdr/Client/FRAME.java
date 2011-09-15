import java.io.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;


public class FRAME extends JFrame {
  JPanel contentPane;
  JPanel Create_Job_Panels = new JPanel();
  JLabel Job_Name_Label = new JLabel();
  JLabel Linux_Executable_Label = new JLabel();
  JLabel Command_Line_Parameters_Label = new JLabel();
  JLabel Frames_Label = new JLabel();
  JLabel User_Name_Label = new JLabel();
  JTextField Job_Name = new JTextField();
  JTextField Linux_Executable = new JTextField();
  JTextField Frames = new JTextField();
  JTextField User_Name = new JTextField();
  JTextField Command_Line_Parameters = new JTextField();
  GridLayout gridLayout1 = new GridLayout();
  JButton Dispatch_Job = new JButton();
  JLabel Server_Name_Label = new JLabel();
  JTextField Server_Name = new JTextField();
  JLabel Server_Directory_Label = new JLabel();
  JTextField Server_Directory = new JTextField();
  BorderLayout borderLayout1 = new BorderLayout();
  JPanel Linux_Executable_Frame = new JPanel();
  JButton Linux_File_Browse = new JButton();
  BorderLayout borderLayout3 = new BorderLayout();
  JPanel Server_Name_Panel = new JPanel();
  JPanel Server_Directory_Panel = new JPanel();
  JPanel Job_Name_Row = new JPanel();
  JPanel jPanel5 = new JPanel();
  JPanel Linux_Executable_Row = new JPanel();
  JPanel Frames_Row = new JPanel();
  JPanel File_Server_Name_Row = new JPanel();
  BorderLayout borderLayout4 = new BorderLayout();
  BorderLayout borderLayout5 = new BorderLayout();
  BorderLayout borderLayout6 = new BorderLayout();
  BorderLayout borderLayout8 = new BorderLayout();
  BorderLayout borderLayout9 = new BorderLayout();
  BorderLayout borderLayout10 = new BorderLayout();
  BorderLayout borderLayout11 = new BorderLayout();
  JTabbedPane jTabbedPane1 = new JTabbedPane();
  JPanel Main_Create_Job_Panel = new JPanel();
  BorderLayout borderLayout12 = new BorderLayout();
  JPanel Add_Slaves_Panel = new JPanel();
  JPanel jPanel10 = new JPanel();
  BorderLayout borderLayout13 = new BorderLayout();
  JPanel jPanel11 = new JPanel();
  GridLayout gridLayout3 = new GridLayout();
  JPanel jPanel12 = new JPanel();
  JPanel jPanel13 = new JPanel();
  GridLayout gridLayout4 = new GridLayout();
  BorderLayout borderLayout14 = new BorderLayout();
  BorderLayout borderLayout15 = new BorderLayout();
  JButton Send_Add_Slaves_Commands = new JButton();
  JPanel jPanel2 = new JPanel();
  JLabel Add_Linux_Slaves_Label = new JLabel();
  JTextField Add_Linux_Slaves = new JTextField();
  BorderLayout borderLayout16 = new BorderLayout();
  JPanel File_Server_Directory_Row = new JPanel();
  JLabel jLabel2 = new JLabel();
  JTextField File_Server_Name = new JTextField();
  JLabel jLabel3 = new JLabel();
  JTextField File_Server_Directory = new JTextField();
  BorderLayout borderLayout17 = new BorderLayout();
  BorderLayout borderLayout18 = new BorderLayout();
  JPanel jPanel14 = new JPanel();
  BorderLayout borderLayout19 = new BorderLayout();
  JPanel jPanel15 = new JPanel();
  BorderLayout borderLayout110 = new BorderLayout();
  JPanel jPanel16 = new JPanel();
  BorderLayout borderLayout111 = new BorderLayout();
  JPanel Command_Line_Row = new JPanel();
  JPanel jPanel18 = new JPanel();
  JButton Add_Add_Spire = new JButton();
  JButton Add_Add_Chromium = new JButton();
  GridLayout gridLayout2 = new GridLayout();
  JLabel jLabel4 = new JLabel();
  JPanel Kill_Frames_Panel = new JPanel();
  JPanel jPanel1 = new JPanel();
  BorderLayout borderLayout2 = new BorderLayout();
  GridLayout gridLayout6 = new GridLayout();
  JPanel jPanel3 = new JPanel();
  BorderLayout borderLayout7 = new BorderLayout();
  JLabel jLabel5 = new JLabel();
  JTextField Kill_Job_ID = new JTextField();
  JPanel jPanel4 = new JPanel();
  JPanel jPanel6 = new JPanel();
  JLabel jLabel6 = new JLabel();
  BorderLayout borderLayout20 = new BorderLayout();
  JTextField Frames_To_Kill = new JTextField();
  JPanel jPanel7 = new JPanel();
  JPanel jPanel8 = new JPanel();
  JPanel jPanel9 = new JPanel();
  JButton Send_Kill_Frames_Command = new JButton();
  JLabel jLabel7 = new JLabel();
  JPanel Retry_Frames = new JPanel();
  BorderLayout borderLayout21 = new BorderLayout();
  JPanel jPanel17 = new JPanel();
  GridLayout gridLayout7 = new GridLayout();
  JPanel jPanel20 = new JPanel();
  JPanel jPanel21 = new JPanel();
  JPanel jPanel22 = new JPanel();
  JPanel jPanel23 = new JPanel();
  JPanel jPanel24 = new JPanel();
  JPanel jPanel25 = new JPanel();
  JButton Send_Retry_Frames_Command = new JButton();
  JLabel jLabel10 = new JLabel();
  BorderLayout borderLayout22 = new BorderLayout();
  BorderLayout borderLayout23 = new BorderLayout();
  BorderLayout borderLayout24 = new BorderLayout();
  BorderLayout borderLayout25 = new BorderLayout();
  BorderLayout borderLayout26 = new BorderLayout();
  BorderLayout borderLayout27 = new BorderLayout();
  JLabel jLabel11 = new JLabel();
  JTextField Retry_Job_ID = new JTextField();
  JTextField Frames_To_Retry = new JTextField();
  JPanel Cancel_Frames = new JPanel();
  JPanel jPanel26 = new JPanel();
  JLabel jLabel12 = new JLabel();
  JPanel jPanel27 = new JPanel();
  JButton Send_Cancel_Frames_Command = new JButton();
  BorderLayout borderLayout28 = new BorderLayout();
  JTextField Frames_To_Cancel = new JTextField();
  JPanel jPanel110 = new JPanel();
  BorderLayout borderLayout29 = new BorderLayout();
  JPanel jPanel28 = new JPanel();
  JTextField Cancel_Job_ID = new JTextField();
  BorderLayout borderLayout210 = new BorderLayout();
  JPanel jPanel29 = new JPanel();
  BorderLayout borderLayout211 = new BorderLayout();
  JPanel jPanel210 = new JPanel();
  BorderLayout borderLayout212 = new BorderLayout();
  BorderLayout borderLayout213 = new BorderLayout();
  GridLayout gridLayout8 = new GridLayout();
  JLabel jLabel13 = new JLabel();
  JPanel jPanel211 = new JPanel();
  GridLayout gridLayout9 = new GridLayout();
  JPanel Remove_Slaves_Panel = new JPanel();
  GridLayout gridLayout10 = new GridLayout();
  JPanel jPanel111 = new JPanel();
  BorderLayout borderLayout112 = new BorderLayout();
  JButton Remove_Add_Spire = new JButton();
  JPanel jPanel112 = new JPanel();
  JLabel jLabel9 = new JLabel();
  JButton Remove_Add_Chromium = new JButton();
  GridLayout gridLayout11 = new GridLayout();
  JPanel jPanel113 = new JPanel();
  JPanel jPanel114 = new JPanel();
  GridLayout gridLayout12 = new GridLayout();
  JPanel jPanel31 = new JPanel();
  BorderLayout borderLayout113 = new BorderLayout();
  BorderLayout borderLayout114 = new BorderLayout();
  JPanel jPanel115 = new JPanel();
  BorderLayout borderLayout115 = new BorderLayout();
  BorderLayout borderLayout116 = new BorderLayout();
  JLabel jLabel14 = new JLabel();
  JTextField Remove_Slaves = new JTextField();
  JPanel jPanel116 = new JPanel();
  JPanel jPanel117 = new JPanel();
  GridLayout gridLayout13 = new GridLayout();
  JLabel Add_Linux_Slaves_Label1 = new JLabel();
  JPanel jPanel118 = new JPanel();
  JButton Send_Remove_Slaves_Commands = new JButton();
  BorderLayout borderLayout117 = new BorderLayout();
  JLabel Slaves_Label = new JLabel();
  JLabel jLabel1 = new JLabel();
  JLabel jLabel8 = new JLabel();
  JTextField Available_Memory = new JTextField();
  JSlider Slave_Priority = new JSlider();
  JPanel jPanel19 = new JPanel();
  JLabel jLabel15 = new JLabel();
  BorderLayout borderLayout30 = new BorderLayout();
  JSlider Job_Priority = new JSlider();
  JPanel jPanel30 = new JPanel();
  JPanel jPanel32 = new JPanel();
  JPanel jPanel33 = new JPanel();
  JPanel jPanel34 = new JPanel();
  JPanel jPanel35 = new JPanel();
  JPanel jPanel36 = new JPanel();
  JLabel jLabel16 = new JLabel();
  BorderLayout borderLayout31 = new BorderLayout();
  JTextField Memory_Requirement = new JTextField();

  //Construct the frame
  public FRAME() {
    enableEvents(AWTEvent.WINDOW_EVENT_MASK);
    try {
      jbInit();
    }
    catch(Exception e) {
      e.printStackTrace();
    }
  }

  //Component initialization
  private void jbInit() throws Exception  {
    contentPane = (JPanel) this.getContentPane();
    contentPane.setLayout(borderLayout1);
    this.setSize(new Dimension(957, 420));
    this.setTitle("PhysBAM Distributed Renderer");
    contentPane.setBackground(SystemColor.control);
    Create_Job_Panels.setLayout(gridLayout1);
    Job_Name_Label.setMaximumSize(new Dimension(200, 15));
    Job_Name_Label.setMinimumSize(new Dimension(200, 15));
    Job_Name_Label.setPreferredSize(new Dimension(200, 15));
    Job_Name_Label.setText("Job Name");
    Linux_Executable_Label.setMaximumSize(new Dimension(200, 15));
    Linux_Executable_Label.setMinimumSize(new Dimension(200, 15));
    Linux_Executable_Label.setPreferredSize(new Dimension(200, 15));
    Linux_Executable_Label.setText("Linux Executable (Local Location)");
    Command_Line_Parameters_Label.setMaximumSize(new Dimension(200, 15));
    Command_Line_Parameters_Label.setMinimumSize(new Dimension(200, 15));
    Command_Line_Parameters_Label.setPreferredSize(new Dimension(200, 15));
    Command_Line_Parameters_Label.setText("Command Line Parameters");
    Frames_Label.setMaximumSize(new Dimension(200, 15));
    Frames_Label.setMinimumSize(new Dimension(200, 15));
    Frames_Label.setPreferredSize(new Dimension(200, 15));
    Frames_Label.setText("Frames (eg: m-n, 3,4,5,8,9)");
    User_Name_Label.setMaximumSize(new Dimension(200, 15));
    User_Name_Label.setMinimumSize(new Dimension(200, 15));
    User_Name_Label.setPreferredSize(new Dimension(200, 15));
    User_Name_Label.setText("User Name");
    Linux_Executable.setMinimumSize(new Dimension(200, 21));
    Linux_Executable.setOpaque(true);
    Linux_Executable.setPreferredSize(new Dimension(200, 21));
    Linux_Executable.setToolTipText("");
    Linux_Executable.setText("");
    gridLayout1.setColumns(1);
    gridLayout1.setHgap(0);
    gridLayout1.setRows(9);
    Dispatch_Job.setFont(new java.awt.Font("Dialog", 1, 11));
    Dispatch_Job.setText("Send Create Job Command");
    Dispatch_Job.addActionListener(new FRAME_Dispatch_Job_actionAdapter(this));
    Server_Name_Label.setMaximumSize(new Dimension(200, 15));
    Server_Name_Label.setMinimumSize(new Dimension(200, 15));
    Server_Name_Label.setPreferredSize(new Dimension(200, 15));
    Server_Name_Label.setText("Server Name");
    Server_Name.setText("joint.stanford.edu");
    Server_Directory_Label.setMaximumSize(new Dimension(200, 15));
    Server_Directory_Label.setMinimumSize(new Dimension(200, 15));
    Server_Directory_Label.setPreferredSize(new Dimension(200, 15));
    Server_Directory_Label.setText("Server Directory");
    Server_Directory.setToolTipText("");
    Server_Directory.setText("/usr/local/pdr");
    Command_Line_Parameters.setText("");
    Frames.setText("");
    User_Name.setText("");
    Linux_File_Browse.setText("Browse");
    Linux_File_Browse.addActionListener(new FRAME_Linux_File_Browse_actionAdapter(this));
    Linux_Executable_Frame.setLayout(borderLayout3);
    Server_Directory_Panel.setLayout(borderLayout4);
    Server_Name_Panel.setLayout(borderLayout5);
    Job_Name_Row.setLayout(borderLayout6);
    jPanel5.setLayout(borderLayout8);
    Linux_Executable_Row.setLayout(borderLayout9);
    Frames_Row.setLayout(borderLayout10);
    File_Server_Name_Row.setLayout(borderLayout11);

    Main_Create_Job_Panel.setDebugGraphicsOptions(0);
    Main_Create_Job_Panel.setLayout(borderLayout12);
    Add_Slaves_Panel.setLayout(borderLayout14);
    jPanel10.setLayout(borderLayout13);
    jPanel11.setLayout(gridLayout3);
    gridLayout3.setColumns(1);
    gridLayout3.setRows(4);
    gridLayout3.setVgap(0);
    jPanel12.setLayout(gridLayout4);
    gridLayout4.setColumns(1);
    gridLayout4.setRows(8);
    jPanel13.setLayout(borderLayout15);
    Send_Add_Slaves_Commands.setFont(new java.awt.Font("Dialog", 1, 11));
    Send_Add_Slaves_Commands.setText("Send Add Slaves Command");
    Send_Add_Slaves_Commands.addActionListener(new FRAME_Send_Add_Slaves_Commands_actionAdapter(this));
    Slaves_Label.setMaximumSize(new Dimension(200, 15));
    Slaves_Label.setMinimumSize(new Dimension(200, 15));
    Slaves_Label.setPreferredSize(new Dimension(200, 15));
    Slaves_Label.setText("Add Linux Slaves");
    jPanel2.setLayout(borderLayout16);
    jLabel2.setMaximumSize(new Dimension(200, 15));
    jLabel2.setMinimumSize(new Dimension(200, 15));
    jLabel2.setPreferredSize(new Dimension(200, 15));
    jLabel2.setText("File Server Name");
    jLabel3.setMaximumSize(new Dimension(200, 15));
    jLabel3.setMinimumSize(new Dimension(200, 15));
    jLabel3.setPreferredSize(new Dimension(200, 15));
    jLabel3.setToolTipText("");
    jLabel3.setText("File Server Directory");
    File_Server_Directory_Row.setLayout(borderLayout17);
    File_Server_Name.setText("");
    File_Server_Directory.setText("");
    jPanel14.setLayout(borderLayout18);
    jPanel15.setLayout(borderLayout19);
    jPanel16.setLayout(borderLayout110);
    Command_Line_Row.setLayout(borderLayout111);
    Add_Add_Spire.setText("Add All Spire");
    Add_Add_Spire.addActionListener(new FRAME_Add_Add_Spire_actionAdapter(this));
    Add_Add_Chromium.setText("Add All Chromium");
    Add_Add_Chromium.addActionListener(new FRAME_Add_Add_Chromium_actionAdapter(this));
    Add_Linux_Slaves.setText("");
    jPanel18.setLayout(gridLayout2);
    jLabel4.setText("        NOTE: If rendering on Chromiums, remember to use an executable " +
    "built for \'pentium\'!");
    Kill_Frames_Panel.setLayout(borderLayout2);
    jPanel1.setLayout(gridLayout6);
    gridLayout6.setColumns(1);
    gridLayout6.setRows(8);
    jPanel3.setMaximumSize(new Dimension(32767, 32767));
    jPanel3.setLayout(borderLayout7);
    jLabel5.setMaximumSize(new Dimension(200, 15));
    jLabel5.setMinimumSize(new Dimension(200, 15));
    jLabel5.setPreferredSize(new Dimension(200, 15));
    jLabel5.setText("Job ID");
    Kill_Job_ID.setText("");
    jLabel6.setMaximumSize(new Dimension(200, 15));
    jLabel6.setMinimumSize(new Dimension(200, 15));
    jLabel6.setPreferredSize(new Dimension(200, 15));
    jLabel6.setText("Frames to Kill");
    jPanel4.setLayout(borderLayout20);
    Frames_To_Kill.setText("");
    Send_Kill_Frames_Command.setFont(new java.awt.Font("Dialog", 1, 11));
    Send_Kill_Frames_Command.setForeground(Color.black);
    Send_Kill_Frames_Command.setText("Send Kill Frames Command");
    Send_Kill_Frames_Command.addActionListener(new FRAME_Send_Kill_Frames_Command_actionAdapter(this));
    jLabel7.setText("                   To kill jobs currently running on a machine, use " +
    "the kill frames command");
    Retry_Frames.setLayout(borderLayout21);
    jPanel17.setLayout(gridLayout7);
    gridLayout7.setColumns(1);
    gridLayout7.setRows(8);
    Send_Retry_Frames_Command.setFont(new java.awt.Font("Dialog", 1, 11));
    Send_Retry_Frames_Command.setText("Send Retry Frames Command");
    Send_Retry_Frames_Command.addActionListener(new FRAME_Send_Retry_Frames_Command_actionAdapter(this));
    jLabel10.setMaximumSize(new Dimension(200, 15));
    jLabel10.setMinimumSize(new Dimension(200, 15));
    jLabel10.setPreferredSize(new Dimension(200, 15));
    jLabel10.setText("Job ID");
    jPanel20.setLayout(borderLayout22);
    jPanel21.setLayout(borderLayout23);
    jPanel22.setLayout(borderLayout24);
    jPanel23.setLayout(borderLayout25);
    jPanel24.setLayout(borderLayout26);
    jPanel25.setLayout(borderLayout27);
    jLabel11.setMaximumSize(new Dimension(200, 15));
    jLabel11.setMinimumSize(new Dimension(200, 15));
    jLabel11.setPreferredSize(new Dimension(200, 15));
    jLabel11.setText("Frames To Retry");
    Frames_To_Retry.setText("");
    Retry_Job_ID.setText("");
    jPanel26.setLayout(borderLayout212);
    jLabel12.setText("Frames To Cancel (Does not kill frames)");
    jLabel12.setPreferredSize(new Dimension(200, 15));
    jLabel12.setMinimumSize(new Dimension(200, 15));
    jLabel12.setMaximumSize(new Dimension(200, 15));
    jPanel27.setLayout(borderLayout213);
    Send_Cancel_Frames_Command.addActionListener(new FRAME_Send_Cancel_Frames_Command_actionAdapter(this));
    Send_Cancel_Frames_Command.setText("Send Cancel Frames Command");
    Send_Cancel_Frames_Command.addActionListener(new FRAME_Send_Cancel_Frames_Command_actionAdapter(this));
    Send_Cancel_Frames_Command.setFont(new java.awt.Font("Dialog", 1, 11));
    Frames_To_Cancel.setText("");
    jPanel110.setLayout(gridLayout8);
    jPanel28.setLayout(borderLayout210);
    Cancel_Job_ID.setText("");
    jPanel29.setLayout(borderLayout211);
    jPanel210.setLayout(borderLayout28);
    gridLayout8.setColumns(1);
    gridLayout8.setRows(8);
    jLabel13.setMaximumSize(new Dimension(200, 15));
    jLabel13.setMinimumSize(new Dimension(200, 15));
    jLabel13.setPreferredSize(new Dimension(200, 15));
    jLabel13.setText("Job ID");
    jPanel211.setLayout(borderLayout29);
    Cancel_Frames.setLayout(gridLayout9);
    Remove_Slaves_Panel.setLayout(gridLayout10);
    jPanel111.setLayout(borderLayout116);
    Remove_Add_Spire.addActionListener(new FRAME_Remove_Add_Spire_actionAdapter(this));
    Remove_Add_Spire.setText("Add All Spire");
    Remove_Add_Spire.addActionListener(new FRAME_Remove_Add_Spire_actionAdapter(this));
    Remove_Add_Spire.setSelectedIcon(null);
    jPanel112.setLayout(gridLayout13);
    jLabel9.setMaximumSize(new Dimension(200, 15));
    jLabel9.setMinimumSize(new Dimension(200, 15));
    jLabel9.setPreferredSize(new Dimension(200, 15));
    jLabel9.setText("Remove Slaves (Doesn\'t End Jobs)");
    Remove_Add_Chromium.setText("Add All Chromium");
    Remove_Add_Chromium.addActionListener(new FRAME_Remove_Add_Chromium_actionAdapter(this));
    Remove_Add_Chromium.addActionListener(new FRAME_Remove_Add_Chromium_actionAdapter(this));
    gridLayout11.setColumns(1);
    gridLayout11.setRows(8);
    jPanel113.setLayout(borderLayout112);
    jPanel114.setLayout(borderLayout113);
    jPanel31.setLayout(borderLayout114);
    jPanel115.setLayout(gridLayout11);
    jLabel14.setText("                   To kill jobs currently running on a machine, use " +
    "the kill frames command");
    Remove_Slaves.setToolTipText("");
    Remove_Slaves.setText("");
    jPanel116.setLayout(gridLayout12);
    jPanel117.setLayout(borderLayout115);
    jPanel118.setLayout(borderLayout117);
    Send_Remove_Slaves_Commands.setFont(new java.awt.Font("Dialog", 1, 11));
    Send_Remove_Slaves_Commands.setText("Send Remove Slaves Command");
    Send_Remove_Slaves_Commands.addActionListener(new FRAME_Send_Remove_Slaves_Commands_actionAdapter(this));
    Send_Remove_Slaves_Commands.addActionListener(new FRAME_Send_Remove_Slaves_Commands_actionAdapter(this));
    Slaves_Label.setMaximumSize(new Dimension(200, 15));
    Slaves_Label.setMinimumSize(new Dimension(200, 15));
    Slaves_Label.setPreferredSize(new Dimension(200, 15));
    Slaves_Label.setText("Add Slaves");
    jLabel1.setMaximumSize(new Dimension(200, 15));
    jLabel1.setMinimumSize(new Dimension(200, 15));
    jLabel1.setPreferredSize(new Dimension(200, 15));
    jLabel1.setText("Available Memory on Slaves (MB)");
    jLabel8.setMaximumSize(new Dimension(200, 15));
    jLabel8.setMinimumSize(new Dimension(200, 15));
    jLabel8.setPreferredSize(new Dimension(200, 15));
    jLabel8.setText("Minimum Job Priority for Slaves (0=min)");
    Available_Memory.setText("");
    Slave_Priority.setOrientation(JSlider.HORIZONTAL);
    Slave_Priority.setMajorTickSpacing(1);
    Slave_Priority.setMaximum(5);
    Slave_Priority.setMinimum(0);
    Slave_Priority.setMinorTickSpacing(1);
    Slave_Priority.setPaintLabels(true);
    Slave_Priority.setPaintTrack(true);
    Slave_Priority.setOpaque(true);
    jLabel15.setMaximumSize(new Dimension(200, 15));
    jLabel15.setMinimumSize(new Dimension(200, 15));
    jLabel15.setPreferredSize(new Dimension(200, 15));
    jLabel15.setText("Job Priority (Starves Lower Priority Jobs)");
    jPanel19.setLayout(borderLayout30);
    Job_Priority.setMajorTickSpacing(1);
    Job_Priority.setMaximum(5);
    Job_Priority.setMinorTickSpacing(0);
    Job_Priority.setPaintLabels(true);
    Job_Priority.setSnapToTicks(true);
    Job_Priority.setValue(2);
    gridLayout9.setRows(1);
    jLabel16.setMaximumSize(new Dimension(200, 15));
    jLabel16.setMinimumSize(new Dimension(200, 15));
    jLabel16.setPreferredSize(new Dimension(200, 15));
    jLabel16.setRequestFocusEnabled(true);
    jLabel16.setText("Memory Requirement (MB)");
    jPanel36.setLayout(borderLayout31);
    Memory_Requirement.setText("");
    Server_Name_Panel.add(Server_Name_Label, BorderLayout.WEST);
    Server_Name_Panel.add(Server_Name, BorderLayout.CENTER);
    Server_Directory_Panel.add(Server_Directory_Label, BorderLayout.WEST);
    Server_Directory_Panel.add(Server_Directory, BorderLayout.CENTER);
    jPanel11.add(jPanel5, null);
    jPanel11.add(Server_Name_Panel, null);
    jPanel11.add(Server_Directory_Panel, null);
    Create_Job_Panels.add(Job_Name_Row, null);
    Job_Name_Row.add(Job_Name_Label, BorderLayout.WEST);
    Job_Name_Row.add(Job_Name, BorderLayout.CENTER);
    jPanel5.add(User_Name_Label, BorderLayout.WEST);
    jPanel5.add(User_Name, BorderLayout.CENTER);
    Create_Job_Panels.add(Linux_Executable_Row, null);
    Linux_Executable_Frame.add(Linux_Executable,  BorderLayout.CENTER);
    Linux_Executable_Frame.add(Linux_File_Browse,  BorderLayout.EAST);
    Create_Job_Panels.add(Frames_Row, null);
    Linux_Executable_Row.add(Linux_Executable_Label, BorderLayout.WEST);
    Command_Line_Row.add(Command_Line_Parameters_Label, BorderLayout.WEST);
    Command_Line_Row.add(Command_Line_Parameters, BorderLayout.CENTER);
    Create_Job_Panels.add(File_Server_Directory_Row, null);
    Create_Job_Panels.add(File_Server_Name_Row, null);
    Create_Job_Panels.add(Command_Line_Row, null);
    Create_Job_Panels.add(jPanel36, null);
    jPanel36.add(jLabel16,  BorderLayout.WEST);
    Create_Job_Panels.add(jPanel19, null);
    jPanel19.add(jLabel15,  BorderLayout.WEST);
    jPanel19.add(Job_Priority, BorderLayout.CENTER);
    Create_Job_Panels.add(Dispatch_Job, null);
    Linux_Executable_Row.add(Linux_Executable_Frame, BorderLayout.CENTER);
    Frames_Row.add(Frames_Label,  BorderLayout.WEST);
    Frames_Row.add(Frames, BorderLayout.CENTER);
    File_Server_Name_Row.add(jLabel2, BorderLayout.WEST);
    File_Server_Name_Row.add(File_Server_Name, BorderLayout.CENTER);
    File_Server_Directory_Row.add(jLabel3, BorderLayout.WEST);
    File_Server_Directory_Row.add(File_Server_Directory, BorderLayout.CENTER);
    contentPane.add(jTabbedPane1,  BorderLayout.CENTER);
    jTabbedPane1.add(Main_Create_Job_Panel,   "Create Job");
    Main_Create_Job_Panel.add(Create_Job_Panels, BorderLayout.CENTER);
    jTabbedPane1.add(Add_Slaves_Panel,    "Add Slaves");
    jTabbedPane1.add(Remove_Slaves_Panel,  "Remove Slaves");
    Remove_Slaves_Panel.add(jPanel115, null);
    jPanel114.add(Add_Linux_Slaves_Label1, BorderLayout.WEST);
    jPanel114.add(jPanel116, BorderLayout.EAST);
    jPanel115.add(jPanel111, null);
    jPanel111.add(jLabel9,  BorderLayout.WEST);
    jPanel111.add(Remove_Slaves,  BorderLayout.CENTER);
    jPanel111.add(jPanel112,  BorderLayout.EAST);
    jPanel115.add(jPanel118, null);
    jPanel115.add(jPanel117, null);
    jPanel115.add(jPanel114, null);
    jPanel112.add(Remove_Add_Chromium, null);
    jPanel112.add(Remove_Add_Spire, null);
    jPanel115.add(jPanel31, null);
    jPanel31.add(jLabel14, BorderLayout.WEST);
    jPanel115.add(jPanel113, null);
    jPanel115.add(jPanel34, null);
    jPanel115.add(Send_Remove_Slaves_Commands, null);
    Add_Slaves_Panel.add(jPanel12,  BorderLayout.CENTER);
    contentPane.add(jPanel11, BorderLayout.NORTH);
    jPanel12.add(jPanel10, null);
    jPanel10.add(Slaves_Label, BorderLayout.WEST);
    jPanel10.add(Add_Linux_Slaves, BorderLayout.CENTER);
    jPanel10.add(jPanel18, BorderLayout.EAST);
    jPanel18.add(Add_Add_Chromium, null);
    jPanel18.add(Add_Add_Spire, null);
    jPanel12.add(jPanel16, null);
    jPanel16.add(jLabel1, BorderLayout.WEST);
    jPanel16.add(Available_Memory, BorderLayout.CENTER);
    jPanel12.add(jPanel15, null);
    jPanel15.add(jLabel8,  BorderLayout.WEST);
    jPanel12.add(jPanel2, null);
    jPanel12.add(jPanel13, null);
    jPanel13.add(jLabel7, BorderLayout.WEST);
    jPanel13.add(jLabel4, BorderLayout.NORTH);
    jTabbedPane1.add(Kill_Frames_Panel,   "Kill Frames");
    Kill_Frames_Panel.add(jPanel1,  BorderLayout.CENTER);
    jPanel1.add(jPanel3, null);
    jPanel3.add(jLabel5,  BorderLayout.WEST);
    jPanel3.add(Kill_Job_ID,  BorderLayout.CENTER);
    jPanel1.add(jPanel4, null);
    jPanel4.add(jLabel6,  BorderLayout.WEST);
    jPanel1.add(jPanel6, null);
    jPanel4.add(Frames_To_Kill, BorderLayout.CENTER);
    jPanel1.add(jPanel7, null);
    jPanel1.add(jPanel8, null);
    jPanel1.add(jPanel9, null);
    jPanel1.add(jPanel33, null);
    jPanel1.add(Send_Kill_Frames_Command, null);
    jTabbedPane1.add(Retry_Frames,   "Retry Frames");
    Retry_Frames.add(jPanel17,  BorderLayout.CENTER);
    jPanel17.add(jPanel20, null);
    jPanel17.add(jPanel21, null);
    jPanel17.add(jPanel22, null);
    jPanel17.add(jPanel23, null);
    jPanel17.add(jPanel24, null);
    jPanel17.add(jPanel25, null);
    jPanel17.add(jPanel32, null);
    jPanel17.add(Send_Retry_Frames_Command, null);
    jTabbedPane1.add(Cancel_Frames,  "Cancel Frames");
    jPanel20.add(jLabel10,  BorderLayout.WEST);
    jPanel21.add(jLabel11, BorderLayout.WEST);
    jPanel20.add(Retry_Job_ID, BorderLayout.CENTER);
    jPanel21.add(Frames_To_Retry, BorderLayout.CENTER);
    jPanel28.add(jLabel13,  BorderLayout.WEST);
    jPanel28.add(Cancel_Job_ID,  BorderLayout.CENTER);
    jPanel110.add(jPanel28, null);
    jPanel110.add(jPanel210, null);
    jPanel210.add(jLabel12,  BorderLayout.WEST);
    jPanel210.add(Frames_To_Cancel,  BorderLayout.CENTER);
    jPanel110.add(jPanel26, null);
    jPanel110.add(jPanel211, null);
    jPanel110.add(jPanel29, null);
    jPanel110.add(jPanel27, null);
    jPanel110.add(jPanel30, null);
    jPanel110.add(Send_Cancel_Frames_Command, null);
    Cancel_Frames.add(jPanel110, null);
    jPanel10.add(Slaves_Label, BorderLayout.WEST);
    jPanel15.add(Slave_Priority, BorderLayout.CENTER);
    jPanel12.add(jPanel14, null);
    jPanel12.add(jPanel35, null);
    jPanel12.add(Send_Add_Slaves_Commands, null);
    jPanel36.add(Memory_Requirement, BorderLayout.CENTER);

    Slave_Priority.setSnapToTicks(true);
    Slave_Priority.setValue(2);
    Restore_All_Fields();
  }

  //File | Exit action performed
  public void jMenuFileExit_actionPerformed(ActionEvent e) {
    Save_All_Fields();
    System.exit(0);
  }


  //Overridden so we can exit when window is closed
  protected void processWindowEvent(WindowEvent e) {
    super.processWindowEvent(e);
    if (e.getID() == WindowEvent.WINDOW_CLOSING) {
     Save_All_Fields();
     System.exit(0);
    }
  }

  void Save_All_Fields()
  {
      try {
        BufferedWriter out = new BufferedWriter(new FileWriter("PDR_Client_Fields"));

        out.write(User_Name.getText()+"\n");
        out.write(Server_Name.getText()+"\n");
        out.write(Server_Directory.getText()+"\n");
        out.write(Job_Name.getText()+"\n");
        out.write(Linux_Executable.getText()+"\n");
        out.write(Frames.getText()+"\n");
        out.write(File_Server_Name.getText()+"\n");
        out.write(File_Server_Directory.getText()+"\n");
        out.write(Command_Line_Parameters.getText()+"\n");

        out.close();
    } catch (IOException e) {
    }
  }

  void Restore_All_Fields()
  {
      try {
        BufferedReader in = new BufferedReader(new FileReader("PDR_Client_Fields"));

        User_Name.setText(in.readLine());
        Server_Name.setText(in.readLine());
        Server_Directory.setText(in.readLine());
        Job_Name.setText(in.readLine());
        Linux_Executable.setText(in.readLine());
        Frames.setText(in.readLine());
        File_Server_Name.setText(in.readLine());
        File_Server_Directory.setText(in.readLine());
        Command_Line_Parameters.setText(in.readLine());

        in.close();
    } catch (IOException e) {
    }
  }

  boolean Is_Number_At_Front_Of_String(String[] str){
    return (str[0].length()>0 && str[0].charAt(0)>='0' && str[0].charAt(0)<='9');
  }

  int Get_First_Number_From_String(String[] str){
    // count the number of digits
    int number_of_digits=0;
    for(int i=0;i<str[0].length();i++){
      if(str[0].charAt(i)>='0' && str[0].charAt(i)<='9') number_of_digits++;
      else break;}
    String first_number=str[0].substring(0,number_of_digits);
    str[0]=str[0].substring(number_of_digits,str[0].length());
    return Integer.parseInt(first_number);
  }

  int Get_Frames_From_String(int[] frames,String string_of_frames)
  {
    int number_of_frames=0;
    boolean last_character_dash=false;
    String frame_string[]=new String[1];
    frame_string[0]=string_of_frames;

    while(Is_Number_At_Front_Of_String(frame_string)){
      // read the number and process it
      int frame_read=Get_First_Number_From_String(frame_string);
      if(last_character_dash){
        int last_index=number_of_frames-1;
        for(int i=frames[last_index]+1;i<=frame_read;i++) frames[number_of_frames++]=i;}
      else frames[number_of_frames++]=frame_read;
      last_character_dash=false;
      // peel off any separator, and process
      while(!Is_Number_At_Front_Of_String(frame_string)){
        if(frame_string[0].length()==0) break;
        if(frame_string[0].charAt(0)=='-'){last_character_dash=true;}
        frame_string[0]=frame_string[0].substring(1,frame_string[0].length());}}
    System.out.print("Frames:  ");
    for(int i=0;i<number_of_frames;i++) System.out.print(frames[i]+"  ");
    System.out.println();
    return number_of_frames;
  }

  void Print_Process_Output_To_Commandline(String string,Process process) throws IOException
  {
    System.out.println("Output From:  \"" + string + "\"");
    BufferedReader proc_in=new BufferedReader(new InputStreamReader(process.getInputStream()));
    String user_input; while((user_input=proc_in.readLine())!=null)System.out.println(user_input);
  }

  void Dispatch_Job_actionPerformed(ActionEvent e) {

    // get all the text fields (any parsing/error checking would be done here)
    String server_name=Server_Name.getText().trim();
    String server_directory=Server_Directory.getText().trim(); // TODO: make sure that there isn't an ending / on the directory name
    String user_name=User_Name.getText().trim();
    String job_name=Job_Name.getText().trim();
    String linux_executable=Linux_Executable.getText().trim();
    String command_line_parameters=Command_Line_Parameters.getText().trim();
    int memory_requirement=Integer.parseInt(Memory_Requirement.getText().trim());
    int priority=Job_Priority.getValue();

    // replace any x: with /cygdrive/x
    linux_executable=linux_executable.replace('\\','/');
    System.out.println(linux_executable);

    // expand out the frames to render
    int frames[]=new int[100000]; // hardcoding this may be bad...
    int number_of_frames=Get_Frames_From_String(frames,Frames.getText());

    // create the job info file (format for the job info file)
    // job name
    // username
    // frames (space separated)
    String job_info_file=(job_name+"\\n"+user_name+"\\n"+priority+"\\n"+memory_requirement+"\\n");
    for(int i=0;i<number_of_frames;i++) job_info_file=job_info_file.concat(frames[i]+" ");
    job_info_file=job_info_file.concat("\\n");

    System.out.println("Job Info File:\n--------------------\n"+job_info_file+"--------------------\n");

    String ssh_string="plink "+user_name+"@"+server_name+" \"umask 0000;cd "+server_directory+"/Commands";
    Process process;
    try{
      // create the appropriate directories on the host
      String create_data_string=ssh_string+"/Create_Data/;mkdir "+user_name+"-"+job_name+"\"";
      process=Runtime.getRuntime().exec(create_data_string);Print_Process_Output_To_Commandline(create_data_string,process);

      String job_info_file_write=ssh_string+"/Create_Data/"+user_name+"-"+job_name+";/bin/echo -ne \\\""+job_info_file+"\\\"> Job_Info\"";
      process=Runtime.getRuntime().exec(job_info_file_write);Print_Process_Output_To_Commandline(job_info_file_write,process);

      if(linux_executable.length()!=0){
        String upload_linux_executable="pscp \""+linux_executable+"\" \""+user_name+"@"+server_name+":"+server_directory+"/Commands/Create_Data/"+user_name+"-"+job_name+"/Linux_Executable\"";
        process=Runtime.getRuntime().exec(upload_linux_executable);Print_Process_Output_To_Commandline(upload_linux_executable,process);}

      String script_header=
          "file_server_input="+File_Server_Name.getText()+":"+File_Server_Directory.getText()+"/Input\\n"+
          "file_server_common="+File_Server_Name.getText()+":"+File_Server_Directory.getText()+"/Common\\n"+
          "file_server_output="+File_Server_Name.getText()+":"+File_Server_Directory.getText()+"/Output\\n"+
          "server_host="+Server_Name.getText()+"\\n"+
          "server_directory="+Server_Directory.getText()+"\\n"+
          "command_line_arguments="+Command_Line_Parameters.getText()+"\\n";
      String upload_script_header=ssh_string+"/Create_Data/"+user_name+"-"+job_name+";/bin/echo -ne \\\""+script_header+"\\\"> Script_Header\"";
      process=Runtime.getRuntime().exec(upload_script_header);Print_Process_Output_To_Commandline(upload_script_header,process);

      // set the permissions on the directory that we just created
      String set_permissions=ssh_string+"/Create_Data/;chmod -R a+rwX "+user_name+"-"+job_name+"\"";
      process=Runtime.getRuntime().exec(set_permissions);Print_Process_Output_To_Commandline(set_permissions,process);

      String touch_create_file=ssh_string+"/Create;touch "+user_name+"-"+job_name+"\"";
      process=Runtime.getRuntime().exec(touch_create_file);Print_Process_Output_To_Commandline(touch_create_file,process);

    }catch(IOException exception) {
      exception.printStackTrace();
      System.out.println("Error when executing commands");}
  }

  boolean Check_Slave_String(String slave_string)
  {
    int at_name=0;
    int at_first_under=1;
    int at_second_under=2;
    int at_processor_count=3;
    int at_space=4;

    if(slave_string.length()==0) return true;

    int current_state=at_name;
    for(int i=0;i<slave_string.length();i++){
      char current_char=slave_string.charAt(i);
      if(current_state==at_name){
        if(Character.isLetter(current_char))            continue;
        if(Character.isDigit(current_char))             continue;
        if(current_char=='-')                           continue;
        if(current_char=='_')                           {current_state=at_first_under;continue;}}
      if(current_state==at_first_under){
        if(current_char=='_')                           {current_state=at_second_under;continue;}}
      if(current_state==at_second_under){
        if(Character.isDigit(current_char))             {current_state=at_processor_count;continue;}}
      if(current_state==at_processor_count){
        if(current_char==' ')                           {current_state=at_space;continue;}}
      if(current_state==at_space){
        if(Character.isLetter(current_char))            {current_state=at_name;continue;}}
      return false;
    }
    if(current_state==at_processor_count) return true;
    return false;
  }

  void Send_Add_Slaves_Commands_actionPerformed(ActionEvent e) {
    // get all the text fields (any parsing/error checking would be done here)
    String server_name=Server_Name.getText().trim();
    String server_directory=Server_Directory.getText().trim(); // TODO: make sure that there isn't an ending / on the directory name
    String user_name=User_Name.getText().trim();
    String add_linux_slaves=Add_Linux_Slaves.getText().trim();
    int available_memory=Integer.parseInt(Available_Memory.getText().trim());
    int slave_priority=Slave_Priority.getValue();
//    String remove_slaves=Remove_Slaves.getText().trim();

    String ssh_string="plink "+user_name+"@"+server_name+" \"umask 0000;cd "+server_directory+"/Commands";
    Process process;
    try{
      // touch the files to remove servers
 /*     if(remove_slaves.length()>0){
        if(!Check_Slave_String(remove_slaves)){System.out.println("Parse Error in Remove_Slaves string. Syntax is: server-name__procs server-name__procs");return;}
        String touch_remove_slaves=ssh_string+"/Remove_Slave"+";touch "+remove_slaves+"\"";
        process=Runtime.getRuntime().exec(touch_remove_slaves);Print_Process_Output_To_Commandline(touch_remove_slaves,process);}
*/
      // parse out the string to add linux servers
      if(!Check_Slave_String(add_linux_slaves)){System.out.println("Parse Error in Linux string. Syntax is: server-name__procs server-name__procs");return;}
      for(int i=0;i<add_linux_slaves.length();i++) if(add_linux_slaves.charAt(i)!=' ' && (i==0 || add_linux_slaves.charAt(i-1)==' ')){
          add_linux_slaves=add_linux_slaves.substring(0,i)+"Linux__"+slave_priority+"__"+available_memory+"__"+add_linux_slaves.substring(i,add_linux_slaves.length());
          i+=7;}

      if(add_linux_slaves.length()>0){
        String touch_add_linux_slaves=ssh_string+"/Add_Slave"+";touch "+add_linux_slaves+"\"";
        process=Runtime.getRuntime().exec(touch_add_linux_slaves);Print_Process_Output_To_Commandline(touch_add_linux_slaves,process);}

    }catch(IOException exception) {
      exception.printStackTrace();
      System.out.println("Error when executing commands");}
  }


  void Linux_File_Browse_actionPerformed(ActionEvent e) {
    JFileChooser file_browser=new JFileChooser();
    if(file_browser.showOpenDialog(this)==JFileChooser.APPROVE_OPTION) Linux_Executable.setText(file_browser.getSelectedFile().toString());
  }

  void Add_Add_Chromium_actionPerformed(ActionEvent e) {
    System.out.println("Adding all the Chromium machines to the list");
    String add_linux_slaves=Add_Linux_Slaves.getText().trim();
    for(int i=1;i<=32;i++) add_linux_slaves=add_linux_slaves+" chromium"+i+"__1"+" chromium"+i+"__2";
    Add_Linux_Slaves.setText(add_linux_slaves);
  }

  void Add_Add_Spire_actionPerformed(ActionEvent e) {
    System.out.println("Adding all the Spire machines to the list");
    String add_linux_slaves=Add_Linux_Slaves.getText().trim();
    for(int i=1;i<=16;i++) add_linux_slaves=add_linux_slaves+" spire-"+i+"__1"+" spire-"+i+"__2";
    Add_Linux_Slaves.setText(add_linux_slaves);
  }

  void Remove_Add_Chromium_actionPerformed(ActionEvent e) {
    System.out.println("Adding all the Chromium machines to the list");
    String remove_slaves=Remove_Slaves.getText().trim();
    for(int i=1;i<=32;i++) remove_slaves=remove_slaves+" chromium"+i+"__1"+" chromium"+i+"__2";
    Remove_Slaves.setText(remove_slaves);
  }

  void Remove_Add_Spire_actionPerformed(ActionEvent e) {
    System.out.println("Adding all the Spire machines to the list");
    String remove_slaves=Remove_Slaves.getText().trim();
    for(int i=1;i<=16;i++) remove_slaves=remove_slaves+" spire-"+i+"__1"+" spire-"+i+"__2";
    Remove_Slaves.setText(remove_slaves);
  }

  void Send_Kill_Frames_Command_actionPerformed(ActionEvent e) {
    String server_name=Server_Name.getText().trim();
    String server_directory=Server_Directory.getText().trim(); // TODO: make sure that there isn't an ending / on the directory name
    String user_name=User_Name.getText().trim();
    String job_id=Kill_Job_ID.getText().trim();
    int frames[]=new int[100000]; // hardcoding this may be bad...
    int number_of_frames=Get_Frames_From_String(frames,Frames_To_Kill.getText());

    String ssh_string="plink "+user_name+"@"+server_name+" \"umask 0000;cd "+server_directory+"/Commands";

    String touch_files="";
    for(int i=0;i<number_of_frames;i++)
    {
      touch_files=touch_files.concat(" "+job_id+"-"+frames[i]);
    }
    Process process;
    try{
      String touch_create_file=ssh_string+"/Kill;touch"+touch_files+"\"";
      process=Runtime.getRuntime().exec(touch_create_file);Print_Process_Output_To_Commandline(touch_create_file,process);
    }catch(IOException exception) {
      exception.printStackTrace();
      System.out.println("Error when executing commands");}
  }

  void Send_Retry_Frames_Command_actionPerformed(ActionEvent e) {
    String server_name=Server_Name.getText().trim();
    String server_directory=Server_Directory.getText().trim(); // TODO: make sure that there isn't an ending / on the directory name
    String user_name=User_Name.getText().trim();
    String job_id=Retry_Job_ID.getText().trim();
    int frames[]=new int[100000]; // hardcoding this may be bad...
    int number_of_frames=Get_Frames_From_String(frames,Frames_To_Retry.getText());

    String ssh_string="plink "+user_name+"@"+server_name+" \"umask 0000;cd "+server_directory+"/Commands";

    String touch_files="";
    for(int i=0;i<number_of_frames;i++)
    {
      touch_files=touch_files.concat(" "+job_id+"-"+frames[i]);
    }
    Process process;
    try{
      String touch_create_file=ssh_string+"/Retry;touch"+touch_files+"\"";
      process=Runtime.getRuntime().exec(touch_create_file);Print_Process_Output_To_Commandline(touch_create_file,process);
    }catch(IOException exception) {
      exception.printStackTrace();
      System.out.println("Error when executing commands");}
  }

  void Send_Cancel_Frames_Command_actionPerformed(ActionEvent e) {
    String server_name=Server_Name.getText().trim();
    String server_directory=Server_Directory.getText().trim(); // TODO: make sure that there isn't an ending / on the directory name
    String user_name=User_Name.getText().trim();
    String job_id=Cancel_Job_ID.getText().trim();
    int frames[]=new int[100000]; // hardcoding this may be bad...
    int number_of_frames=Get_Frames_From_String(frames,Frames_To_Cancel.getText());

    String ssh_string="plink "+user_name+"@"+server_name+" \"umask 0000;cd "+server_directory+"/Commands";

    String touch_files="";
    for(int i=0;i<number_of_frames;i++)
    {
      touch_files=touch_files.concat(" "+job_id+"-"+frames[i]);
    }
    Process process;
    try{
      String touch_create_file=ssh_string+"/Cancel;touch"+touch_files+"\"";
      process=Runtime.getRuntime().exec(touch_create_file);Print_Process_Output_To_Commandline(touch_create_file,process);
    }catch(IOException exception) {
      exception.printStackTrace();
      System.out.println("Error when executing commands");}
  }

  void Send_Remove_Slaves_Commands_actionPerformed(ActionEvent e) {
   // get all the text fields (any parsing/error checking would be done here)
    String server_name=Server_Name.getText().trim();
    String server_directory=Server_Directory.getText().trim(); // TODO: make sure that there isn't an ending / on the directory name
    String user_name=User_Name.getText().trim();
    String remove_slaves=Remove_Slaves.getText().trim();

    String ssh_string="plink "+user_name+"@"+server_name+" \"umask 0000;cd "+server_directory+"/Commands";
    Process process;
    try{
      // touch the files to remove servers
      if(remove_slaves.length()>0){
        if(!Check_Slave_String(remove_slaves)){System.out.println("Parse Error in Remove_Slaves string. Syntax is: server-name__procs server-name__procs");return;}
        String touch_remove_slaves=ssh_string+"/Remove_Slave"+";touch "+remove_slaves+"\"";
        process=Runtime.getRuntime().exec(touch_remove_slaves);Print_Process_Output_To_Commandline(touch_remove_slaves,process);}

    }catch(IOException exception) {
      exception.printStackTrace();
      System.out.println("Error when executing commands");}
  }
}

class FRAME_Dispatch_Job_actionAdapter implements java.awt.event.ActionListener {
  FRAME adaptee;

  FRAME_Dispatch_Job_actionAdapter(FRAME adaptee) {
    this.adaptee = adaptee;
  }
  public void actionPerformed(ActionEvent e) {
    adaptee.Dispatch_Job_actionPerformed(e);
  }
}

class FRAME_Linux_File_Browse_actionAdapter implements java.awt.event.ActionListener {
  FRAME adaptee;

  FRAME_Linux_File_Browse_actionAdapter(FRAME adaptee) {
    this.adaptee = adaptee;
  }
  public void actionPerformed(ActionEvent e) {
    adaptee.Linux_File_Browse_actionPerformed(e);
  }
}

class FRAME_Send_Add_Slaves_Commands_actionAdapter implements java.awt.event.ActionListener {
  FRAME adaptee;

  FRAME_Send_Add_Slaves_Commands_actionAdapter(FRAME adaptee) {
    this.adaptee = adaptee;
  }
  public void actionPerformed(ActionEvent e) {
    adaptee.Send_Add_Slaves_Commands_actionPerformed(e);
  }
}

class FRAME_Add_Add_Chromium_actionAdapter implements java.awt.event.ActionListener {
  FRAME adaptee;

  FRAME_Add_Add_Chromium_actionAdapter(FRAME adaptee) {
    this.adaptee = adaptee;
  }
  public void actionPerformed(ActionEvent e) {
    adaptee.Add_Add_Chromium_actionPerformed(e);
  }
}

class FRAME_Add_Add_Spire_actionAdapter implements java.awt.event.ActionListener {
  FRAME adaptee;

  FRAME_Add_Add_Spire_actionAdapter(FRAME adaptee) {
    this.adaptee = adaptee;
  }
  public void actionPerformed(ActionEvent e) {
    adaptee.Add_Add_Spire_actionPerformed(e);
  }
}

class FRAME_Send_Kill_Frames_Command_actionAdapter implements java.awt.event.ActionListener {
  FRAME adaptee;

  FRAME_Send_Kill_Frames_Command_actionAdapter(FRAME adaptee) {
    this.adaptee = adaptee;
  }
  public void actionPerformed(ActionEvent e) {
    adaptee.Send_Kill_Frames_Command_actionPerformed(e);
  }
}

class FRAME_Send_Retry_Frames_Command_actionAdapter implements java.awt.event.ActionListener {
  FRAME adaptee;

  FRAME_Send_Retry_Frames_Command_actionAdapter(FRAME adaptee) {
    this.adaptee = adaptee;
  }
  public void actionPerformed(ActionEvent e) {
    adaptee.Send_Retry_Frames_Command_actionPerformed(e);
  }
}

class FRAME_Send_Cancel_Frames_Command_actionAdapter implements java.awt.event.ActionListener {
  FRAME adaptee;

  FRAME_Send_Cancel_Frames_Command_actionAdapter(FRAME adaptee) {
    this.adaptee = adaptee;
  }
  public void actionPerformed(ActionEvent e) {
    adaptee.Send_Cancel_Frames_Command_actionPerformed(e);
  }
}

class FRAME_Send_Remove_Slaves_Commands_actionAdapter implements java.awt.event.ActionListener {
  FRAME adaptee;

  FRAME_Send_Remove_Slaves_Commands_actionAdapter(FRAME adaptee) {
    this.adaptee = adaptee;
  }
  public void actionPerformed(ActionEvent e) {
    adaptee.Send_Remove_Slaves_Commands_actionPerformed(e);
  }
}

class FRAME_Remove_Add_Spire_actionAdapter implements java.awt.event.ActionListener {
  FRAME adaptee;

  FRAME_Remove_Add_Spire_actionAdapter(FRAME adaptee) {
    this.adaptee = adaptee;
  }
  public void actionPerformed(ActionEvent e) {
    adaptee.Remove_Add_Spire_actionPerformed(e);
  }
}

class FRAME_Remove_Add_Chromium_actionAdapter implements java.awt.event.ActionListener {
  FRAME adaptee;

  FRAME_Remove_Add_Chromium_actionAdapter(FRAME adaptee) {
    this.adaptee = adaptee;
  }
  public void actionPerformed(ActionEvent e) {
    adaptee.Remove_Add_Chromium_actionPerformed(e);
  }
}

