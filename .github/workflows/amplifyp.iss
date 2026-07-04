[Setup]
AppName=AmplifyP
AppVersion={#Version}
AppPublisher=Fufu Fang
AppPublisherURL=https://github.com/fangfufu/AmplifyP
AppSupportURL=https://github.com/fangfufu/AmplifyP/issues
AppUpdatesURL=https://github.com/fangfufu/AmplifyP/releases
DefaultDirName={localappdata}\Programs\AmplifyP
DefaultGroupName=AmplifyP
SetupIconFile=src/assets/icon.ico
UninstallDisplayIcon={app}\AmplifyP.exe
LicenseFile=LICENSE
SourceDir=..\..
OutputBaseFilename=amplifyp-windows-setup-{`#Version`}
Compression=lzma2/ultra64
SolidCompression=yes
WizardStyle=modern
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog

[Files]
Source: "build\AmplifyP\*"; DestDir: "{app}"; Flags: recursesubdirs

[Icons]
Name: "{group}\AmplifyP"; Filename: "{app}\AmplifyP.exe"
Name: "{autodesktop}\AmplifyP"; Filename: "{app}\AmplifyP.exe"

[Run]
Filename: "{app}\AmplifyP.exe"; Description: "Launch AmplifyP"; Flags: nowait postinstall skipifsilent
